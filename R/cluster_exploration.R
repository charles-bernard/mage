# mage R package

# a series of functions wrapped to explore each cluster
# as a tree

# recquired internal function 1:
# compute the centroid of the population and for each cluster,
# define the point which is the nearest to this centroid
# ------------------------------------------------------------------
get_starting_pt <- function(scores_mat, partition) {
  # This is the centroid of the population
  centroid_coord <- colMeans(scores_mat);

  clusters <- sort(unique(partition));
  clust_starting_pt <- list();
  for(i in 1:length(clusters)) {
    # Find closest point to centroid in the current cluster
    # This will further be used as the starting point for the cluster exploration
    # ---------------------------------------------------------------
    alldist <- as.matrix(dist(rbind(centroid_coord, scores_mat[which(partition == clusters[i]), ])));
    mindist_to_centroid <- sort(alldist[1,-1], index.return = T)$ix;
    clust_starting_pt[[i]] <- mindist_to_centroid[1];
  }
  names(clust_starting_pt) <- clusters;
  return(clust_starting_pt);
}

# recquired internal function 2:
# compute the ElPiGraph Tree of each cluster to capture
# the different trajectories of its points in the scores space
# ------------------------------------------------------------------
construct_tree <- function(scores_mat, partition) {
  clusters <- sort(unique(partition));

  clusters_tree <- list();
  for(i in 1:length(clusters)) {
    curr_clust <- clusters[i]
    curr_clust_data <- scores_mat[which(partition == curr_clust), ];

    # Build a list of ElPiGraph conservative trees with bootstrapping
    # ---------------------------------------------------------------
    conservative_trees <- computeElasticPrincipalTree(
      X = curr_clust_data,
      # Do not proceed to dimensionality reduction
      Do_PCA = FALSE,
      Mu = 1, Lambda = .01,
      NumNodes = 30, MaxNumberOfIterations = 50,
      # Bootstrapping Parameters
      # nReps = 10, ProbPoint = .9,
      # Avoid default graphical output}
      drawAccuracyComplexity = FALSE, drawEnergy = FALSE, drawPCAView = FALSE);

    # Last tree of the list is the average of all bootstrapped trees
    # ----------------------------------------------------------------
    clusters_tree[[i]] <- conservative_trees[[length(conservative_trees)]];
  }

  names(clusters_tree) <- clusters;
  return(clusters_tree);
}

# recquired internal function 3:
# get info about the trajectories and the nodes of each tree
# and most importantly about the points these objects are associated to
# ------------------------------------------------------------------
get_tree_info <- function(scores_mat, partition, trees, trees_starting_pt) {
  clusters <- sort(unique(partition));
  clust_tree_info <- list();

  for(i in 1:length(clusters)) {

    curr_clust <- clusters[i]
    curr_clust_data <- scores_mat[which(partition == curr_clust), ];
    curr_tree <- trees[[i]];

    clust_tree_info[[i]] <- list();

    # Associate points to nodes
    # --------------------------------------------------------------
    pt_nodes <-
      PartitionData(X = curr_clust_data, NodePositions = curr_tree$NodePositions)$Partition;

    # Identify root node
    # --------------------------------------------------------------
    root <- pt_nodes[trees_starting_pt[[i]]];

    # igraph network from the ElPiGraph tree structure
    # --------------------------------------------------------------
    tree_graph <-
      ConstructGraph(PrintGraph = curr_tree);

    # Get all trajectories of the tree
    # --------------------------------------------------------------
    tree_all_traj <-
      GetSubGraph(Net = tree_graph, Structure = 'end2end');

    # Get trajectories with root node within
    # --------------------------------------------------------------
    tree_traj <- tree_all_traj[sapply(tree_all_traj,
                                      function(x) {any(x[c(1:length(x))] == root)})];

    # If root node is not a leaf, find nearest leaf from it:
    # This will be set as the new root
    # --------------------------------------------------------------
    root <- tree_traj[[which.min(unlist(lapply(tree_traj,
                                               function(x) { which(x == root) })))]][1];


    # If trajectory ends with the root node, revert its node order
    # --------------------------------------------------------------
    tree_traj <- lapply(tree_traj, function(x) {
      if(x[1] == root) { return(x) } else { return(rev(x)) }});
    names(tree_traj) <- 1:length(tree_traj);

    # Associate nodes to all trajectories they are on
    # --------------------------------------------------------------
    node_traj <- list();
    for(n in 1:vcount(tree_graph)) {
      trajs <- NULL;
      for(t in 1:length(tree_traj)) {
        if(n %in% tree_traj[[t]]) {
          trajs <- c(trajs, t);
        }
      }
      node_traj[[n]] <- trajs;
    }

    clust_tree_info[[i]]$root <- root;
    clust_tree_info[[i]]$tree_graph <- tree_graph;
    clust_tree_info[[i]]$tree_traj <- tree_traj;
    clust_tree_info[[i]]$node_traj <- node_traj;
    clust_tree_info[[i]]$pt_nodes <- pt_nodes;
  }

  names(clust_tree_info) <- clusters;
  return(clust_tree_info);
}

# recquired internal function 4:
# Create a tabular file, which will serve as a template to help
# the user make a comprehensive classification of all the significant
# associations in his/her dataset
# ------------------------------------------------------------------
create_classification_template <- function(outdir, trees_info) {
  n_clusters <- length(trees_info);
  clusters <- names(trees_info);
  clusters_vector <- branches_vector <- nodes_vector <- NULL;

  for(i in 1:n_clusters) {
    n_nodes <- length(V(trees_info[[i]]$tree_graph));
    clusters_vector <- c(clusters_vector, rep(clusters[i], n_nodes));
    branches_vector <- c(branches_vector, trees_info[[i]]$node_traj);
    nodes_vector <- c(nodes_vector, as.character(1:n_nodes));
  }
  class_vector <- rep("", length(nodes_vector));

  class_template <- data.table(
    `Cluster` = clusters_vector,
    `Branch` = branches_vector,
    `Node` = nodes_vector,
    `Classification` = class_vector);

  fwrite(class_template,
         file = file.path(outdir, 'associations_classification_template.csv'),
         sep = ",", col.names = TRUE);

  sprintf("A classification file has been produced at %s",
          file.path(outdir, 'associations_classification_template.csv'));
}

# recquired internal function 5:
# A simple function for the creation of directories
# returns subdirectory path
# ------------------------------------------------------------------
create_subdir <- function(main_dir, subdirname) {
  subdir = file.path(main_dir, subdirname);
  if(file.exists(subdir))
    unlink(subdir, recursive = TRUE);
  dir.create(subdir);
  return(subdir);
}

# recquired internal function 6:
# This is the core function
# Visit the tree of the current cluster from root trajectory to child trajectories
# Starting from the node which is the nearest to the centroid of the population
# At each visited node, the function creates a directory
# in which scatterplots of associated points are produced
# It results from this function an arboresence of directories which
# match the trajectory ramifications of the ElPiGraph Tree
# ------------------------------------------------------------------
create_and_fill_traj <-
  function(scores_tab, scores_mat, gene_exp_mat,
           node = node,
           parent_node = -1, parent_trajs = -1, parent_dir,
           tree_info, tree, sampling_ix) {

  # Get all the trajectories of the current node
  # -----------------------------------------------------------------
  node_trajs <- tree_info$node_traj[[node]];

  # Get indexes of the points projected on the current node
  # -----------------------------------------------------------------
  node_ix <- which(tree_info$pt_nodes == node);

  # if the current node doesn't follow all the same trajectories
  # as the parent node (-> we are after a bifurcation)
  # -----------------------------------------------------------------
  if(!identical(node_trajs, parent_trajs)) {
    # get name of the new trajectory(ies)
    # ---------------------------------------------------------------
    traj_name <- paste("traj", paste(node_trajs, collapse = "_"), sep = "_");
    parent_dir <- create_subdir(parent_dir, traj_name);

    # Highlight the new trajectory(ies) on the tree
    # -----------------------------------------------------------------
    pdf(file.path(parent_dir, paste("tree_", traj_name, ".pdf", sep = "")));

    points_on_traj <- unlist(lapply(
      tree_info$node_traj[tree_info$pt_nodes], function(x) any(x %in% node_trajs)));
    pp1 <- PlotPG(X = scores_mat,
                  Main = paste("Points in", traj_name, "\n"),
                  TargetPG = tree,
                  GroupsLab = points_on_traj,
                  NodeLabels = 1:nrow(tree$NodePositions),
                  LabMult = 4, DimToPlot = 1:2,
                  PlotProjections = "onNodes", p.alpha = .5);

    nodes_on_traj <- unlist(lapply(
      tree_info$node_traj, function(x) any(x %in% node_trajs)));
    pp2 <- PlotPG(X = scores_mat,
                  Main = paste(traj_name, "\n"),
                  TargetPG = tree,
                  PGCol = nodes_on_traj,
                  NodeLabels = 1:nrow(tree$NodePositions),
                  LabMult = 4, DimToPlot = 1:2,
                  PlotProjections = "onNodes", p.alpha = .5);

    plot(pp1[[1]]); plot(pp2[[1]]);
    junk <- dev.off();
  }

  if(node %in% names(sampling_ix)) {
    # Get position of the node relative to its trajectories
    # (The node should be at the same position whatever the trajectory it is on)
    # -----------------------------------------------------------------
    node_position_along_traj <- which(tree_info$tree_traj[[node_trajs[1]]] == node);
    dirname <- paste("Pos", node_position_along_traj, "_Node", node, sep = "");
    node_dir <- create_subdir(parent_dir, dirname);

    # scatterplot of sampled individuals within the node
    # -----------------------------------------------------------------
    sampled_node_ix <- sampling_ix[[which(names(sampling_ix) == node)]][, 1];
    if(length(sampled_node_ix) > 0) {
      if(any(sampled_node_ix, na.rm = TRUE)) {
        scatterplot(node_dir, node_ix[sampled_node_ix], scores_tab, gene_exp_mat);
      }
    }

    # Highlight points projected on the node
    # -----------------------------------------------------------------
    pdf(file.path(node_dir, paste("0_", dirname, ".pdf", sep = "")));
    pp3 <- PlotPG(X = scores_mat,
                  Main = dirname,
                  TargetPG = tree,
                  GroupsLab = tree_info$pt_nodes %in% node,
                  NodeLabels = 1:nrow(tree$NodePositions),
                  LabMult = 4, DimToPlot = 1:2,
                  PlotProjections = "onNodes", p.alpha = .5);
    plot(pp3[[1]]);
    junk <- dev.off();
  }

  # Get adjacent node(s) to current node:
  # -----------------------------------------------------------------
  adj_nodes <- adjacent_vertices(tree_info$tree_graph, node)[[1]];

  # Recursivity on adjacent nodes
  # -----------------------------------------------------------------
  # Get rid of already visited parent in the adj_nodes
  if(parent_node != -1) {
    adj_nodes <- adj_nodes[-which(adj_nodes == parent_node)];
  }
  # Do nothing more if current node is a leaf (and not a root)
  # -----------------------------------------------------------------
  if(length(adj_nodes) > 0) {
    for(n in 1:length(adj_nodes)) {
      create_and_fill_traj(
        scores_tab, scores_mat, gene_exp_mat,
        node = adj_nodes[n],
        parent_node = node, parent_trajs = node_trajs, parent_dir = parent_dir,
        tree_info, tree, sampling_ix);
    }
  }
}

# Wrapper of all the previous internal functions
# -----------------------------------------------------------------
#' @title explore_clusters
#'
#' @description Given a table of scores and its partition into k clusters,
#' the function \code{explore_clusters} proposes to explore the types
#' of associations lying within each cluster.
#'
#' To do so, the function measures the dispersion of the associations of each
#' cluster in the space of the scores and depicts this heterogeneity
#' as a tree structure (using ElPiGraph.R package).
#' For instance, nodes along a same branch will be expected to exhibit
#' similar types of associations.
#'
#' Then, the function will visit the tree of each cluster to create an arborescence of
#' directories starting from the root (the root is defined
#' as the branch which is the closest to the centroid of the table of scores)
#' and expanding from adjacent branches to adjacent branches.
#'
#' Inside each directory corresponding to a branch (or a branching point),
#' the function will produce the scatterplots of a few sampled assocations in order
#' for the user to get a glimpse of what type of associations each branch is related to.
#'
#' By visiting this arborescence of directories and visualizing
#' some of the scatterplots of the associations each branch is related to,
#' it enables one to somehow "wander" along a cluster's tree and appreciate the
#' "evolution" of the "shape/type" of the associations from branches to branches
#'
#' @param scores_tab data.table; the table of scores
#' @param partition categorial vector; the partition of the \code{scores_tab} into k clusters.
#' \code{length(partition)} must be equal to \code{nrow(scores_tab)}
#' @param gene_exp_mat matrix: the gene expression matrix used for computing the scores
#' (rownames corresponding to genes and colnames to cells)
#' @param sampling_method method to be used for selecting a few associations in each branch.
#'
#' Methods available:
#' \itemize{
#'   \item{\code{"downsampling"}}
#'   \item{\code{"lowsampling"}}
#'   \item{\code{"upsampling"}}
#'   \item{\code{"centroids"}}
#'   \item{\code{"medoids"}}
#'  }
#'
#' Method recommanded: \code{"downsampling"}
#' @param target_ratio number; fraction of associations to be sampled in each tree's node
#' . Argument \code{target_ratio} will be ignored if argument \code{target_nb} is provided
#' @param target_nb number of associations to be sampled in each tree's node
#' @param output_directory where to create the arborescence of directories for each cluster
#'
#' @importFrom ElPiGraph.R computeElasticPrincipalTree PlotPG ConstructGraph GetSubGraph PartitionData
#' @importFrom igraph adjacent_vertices vcount V
#' @importFrom stats dist
#' @importFrom data.table data.table fwrite
explore_clusters <- function(scores_tab, partition, gene_exp_mat,
                             sampling_method = "downsampling",
                             target_ratio = .3, target_nb = NULL,
                             output_directory = ".") {

  # Create a matrix of scores (SPEARMAN rho is ignored) from the table of scores
  # -----------------------------------------------------------------
  scores_mat <- data.matrix(scores_tab[, 3:(ncol(scores_tab)-1)]);

  # Identify clusters from partition
  # -----------------------------------------------------------------
  clusters <- sort(unique(partition));

  # Build a robust tree (resulting averaged tree from bootstrapping)
  # that captures the heterogeneity of each cluster
  # -----------------------------------------------------------------
  trees <- construct_tree(scores_mat, partition);

  # Get the root point of each tree
  # -----------------------------------------------------------------
  trees_starting_pt <- get_starting_pt(scores_mat, partition);

  # Get relationships between trajectories, nodes and points of each tree
  # -----------------------------------------------------------------
  trees_info <- get_tree_info(scores_mat, partition, trees,  trees_starting_pt);

  # Sample a few points within each node of each tree
  # -----------------------------------------------------------------
  trees_sampling <- list();
  for(i in 1:length(clusters)) {
    trees_sampling[[i]] <-
      sample(scores_tab[partition == clusters[i], ],
             partition = trees_info[[i]]$pt_nodes,
             methods = sampling_method,
             target_ratio = target_ratio, target_nb = target_nb);
  }

  # This is where the tree exploration by arborescence begins
  # -----------------------------------------------------------------
  for(i in 1:length(clusters)) {

    # 1. Create a directory for the current cluster
    # ---------------------------------------------------------------
    curr_cluster_dir <- create_subdir(output_directory, paste("Cluster", clusters[i]));

    # 2. Find the node whose starting point of the tree is associated to
    # ---------------------------------------------------------------
    curr_root_node <- trees_info[[i]]$pt_nodes[trees_starting_pt[[i]]];

    # 3. Plot the tree of the cluster at the root of the cluster directory
    # ---------------------------------------------------------------
    pdf(file.path(curr_cluster_dir, "cluster_tree.pdf"));
    p <- PlotPG(X = scores_mat[partition == clusters[i], ], TargetPG = trees[[i]],
                GroupsLab = trees_info[[i]]$pt_br_brpt,
                NodeLabels = 1:nrow(trees[[i]]$NodePositions),
                LabMult = 4, DimToPlot = 1:2,
                PlotProjections = "onNodes", p.alpha = .5);
    plot(p[[1]]);
    junk <- dev.off();

    # Create the template csv file to help classifying the associations
    # -----------------------------------------------------------------
    create_classification_template(output_directory, trees_info);

    # 4. Create arboresence of directories for the trajectories,
    # meanwhile plotting sampled associations in each node directory
    # (core function)
    # ---------------------------------------------------------------
    create_and_fill_traj (
      scores_tab[partition == clusters[i], ],
      scores_mat[partition == clusters[i], ],
      gene_exp_mat,
      node = curr_root_node,
      parent_node = -1, parent_trajs= -1,
      parent_dir = curr_cluster_dir,
      tree_info = trees_info[[i]],
      tree = trees[[i]],
      sampling_ix = trees_sampling[[i]]);
  }

  return(trees_info);
}
#
# output_directory = "~/Documents/vignettes_out/";
# partition = partition;
# gene_exp_mat = gene_exp_mat;
# sampling_method = "downsampling";
# target_ratio = .8;
# scores_tab = signif_scores;

