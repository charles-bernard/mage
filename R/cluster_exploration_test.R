# mage R package

# a series of functions wrapped to explore each cluster
# as a tree

# recquired internal function 1:
# compute the centroid of the population and for each cluster,
# define the point which is the nearest to this centroid
# ------------------------------------------------------------------
get_clust_starting_pt <- function(scores_mat, partition) {
  # This is the centroid of the population
  centroid_coord <- colMeans(scores_mat);

  clusters <- sort(unique(partition));
  clusters_starting_pt <- list();
  for(i in 1:length(clusters)) {
    # Find closest point to centroid in the current cluster
    # This will further be used as the starting point for the cluster exploration
    alldist <- as.matrix(dist(rbind(centroid_coord, scores_mat[which(partition == clusters[i]), ])));
    mindist_to_centroid <- sort(alldist[1,-1], index.return = T)$ix;
    clusters_starting_pt[[i]] <- mindist_to_centroid[1];
  }
  names(clusters_starting_pt) <- clusters;
  return(clusters_starting_pt)
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
    conservative_trees <- computeElasticPrincipalTree(
      X = curr_clust_data,
      # Do not proceed to dimensionality reduction
      Do_PCA = FALSE,
      Mu = 1, Lambda = .01,
      NumNodes = 30, MaxNumberOfIterations = 50,
      # Bootstrapping Parameters
      nReps = 10, ProbPoint = .9,
      # Avoid default graphical output}
      drawAccuracyComplexity = FALSE, drawEnergy = FALSE, drawPCAView = FALSE);

    # Last tree of the list is the average of all bootstrapped trees
    clusters_tree[[i]] <- conservative_trees[[length(conservative_trees)]];
  }

  names(clusters_tree) <- clusters;
  return(clusters_tree);
}

# recquired internal function 3:
# get info about the branches, the branching points and the nodes of each tree
# and most importantly about the points these structures are associated to
# ------------------------------------------------------------------
get_tree_info <- function(scores_mat, partition, trees) {
  clusters <- sort(unique(partition));
  clusters_tree_info <- list();

  for(i in 1:length(clusters)) {
    curr_clust <- clusters[i]
    curr_clust_data <- scores_mat[which(partition == curr_clust), ];
    curr_tree <- trees[[i]];

    clusters_tree_info[[i]] <- list();
    # igraph network from the ElPiGraph tree structure
    tree_graph <-
      ConstructGraph(PrintGraph = curr_tree);
    # Get branches and branching points of the tree
    tree_br_brpt <-
      GetSubGraph(Net = tree_graph, Structure = 'branches&bpoints');
    # Associate nodes to branches or branching points
    node_br_brpt <- rep("", vcount(tree_graph));
    for(j in 1:length(tree_br_brpt)) {
      node_br_brpt[tree_br_brpt[[j]]] <- names(tree_br_brpt)[j];
    }
    # Associate points to nodes
    pt_nodes <-
      PartitionData(X = curr_clust_data, NodePositions = curr_tree$NodePositions)$Partition;
    # Associate points to branches or branching points
    pt_br_brpt <- node_br_brpt[pt_nodes];

    clusters_tree_info[[i]]$tree_graph <- tree_graph; rm(tree_graph);
    clusters_tree_info[[i]]$tree_br_brpt <- tree_br_brpt; rm(tree_br_brpt);
    clusters_tree_info[[i]]$node_br_brpt <- node_br_brpt; rm(node_br_brpt);
    clusters_tree_info[[i]]$pt_nodes <- pt_nodes; rm(pt_nodes);
    clusters_tree_info[[i]]$pt_br_brpt <- pt_br_brpt; rm(pt_br_brpt);
  }

  names(clusters_tree_info) <- clusters;
  return(clusters_tree_info);
}

get_tree_trajectories <- function(scores_mat, partition,
                                  trees, trees_starting_nodes,
                                  trees_info) {
  clusters <- sort(unique(partition));
  clusters_traj_info <- list();

  for(i in 1:length(clusters)) {
    clusters_traj_info[[i]] <- list();
    tree_e2e <- GetSubGraph(Net = trees_info[[i]]$tree_graph, Structure = 'end2end');
    root <- trees_starting_nodes[[i]];
    selpaths <- tree_e2e[sapply(tree_e2e, function(x){any(x[c(1, length(x))] == root)})];
    selpaths <- lapply(selpaths, function(x) {
      if(x[1] == root){
        return(x)
      } else {
        return(rev(x))
      }
    })
    clusters_traj_info[[i]]$tree_e2e <- tree_e2e;
    clusters_traj_info[[i]]$paths <- selpaths;
  }
  return(clusters_traj_info);
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
    branches_vector <- c(branches_vector, trees_info[[i]]$node_br_brpt);
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
}

# recquired internal function 5:
# A simple function for the creation of directories
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
# Visit the tree of the current cluster from branches to adjacent branches
# Starting from the node which is the nearest to the centroid of the population
# At each branch or branching points, the function creates a directory
# as well as node subdirectories in which scatterplots of associated points are produced
# It results from this function an arboresence of directories
# identical to the one of the ElPiGraph Tree
# ------------------------------------------------------------------
create_and_fill_arborescence <- function(scores_tab, scores_mat, gene_exp_mat,
                                         node = node,
                                         parent_node = -1, parent_struct = -1, parent_dir,
                                         tree_info, tree, sampling_ix) {
  # Get branch or branching point of the current node
  node_struct <- tree_info$node_br_brpt[node];

  # Get indexes of the points projected on the current node
  node_ix <- which(tree_info$pt_nodes == node);

  # if we are visiting a different branch from parent's,
  # create directory with name of current branch
  if(node_struct != parent_struct) {
    parent_dir <- create_subdir(parent_dir, node_struct);

    pdf(file.path(parent_dir, paste("tree_", node_struct, ".pdf", sep = "")));
    pp1 <- PlotPG(X = scores_mat,
                  TargetPG = tree,
                  GroupsLab = tree_info$pt_br_brpt %in% node_struct,
                  NodeLabels = 1:nrow(tree$NodePositions),
                  LabMult = 4, DimToPlot = 1:2, PlotProjections = "onNodes", p.alpha = .5);
    pp2 <- PlotPG(X = scores_mat,
                  TargetPG = tree,
                  PGCol = V(tree_info$tree_graph) %in%
                    tree_info$tree_br_brpt[[which(names(tree_info$tree_br_brpt) == node_struct)]],
                  NodeLabels = 1:nrow(tree$NodePositions),
                  LabMult = 4, DimToPlot = 1:2, PlotProjections = "onNodes", p.alpha = .5);
    plot(pp1[[1]]); plot(pp2[[1]]);
    junk <- dev.off();
  }

  if(node %in% names(sampling_ix)) {
    # Get position of the node relative to its branch (create dir accordingly)
    node_position_along_struct <- which(tree_info$tree_br_brpt[[node_struct]] == node);
    dirname <- paste("Pos", node_position_along_struct, "_Node", node, sep = "");
    node_dir <- create_subdir(parent_dir, dirname);

    # scatterplot of sampled individuals within the node
    sampled_node_ix <- sampling_ix[[which(names(sampling_ix) == node)]][, 1];
    if(any(sampled_node_ix)) {
      scatterplot(node_dir, node_ix[sampled_node_ix], scores_tab, gene_exp_mat);
    }
  }

  # Get adjacent node(s) to current node:
  adj_nodes <- adjacent_vertices(tree_info$tree_graph, node)[[1]];

  # Recursivity on adjacent nodes
  # -----------------------------------------------------------------
  # Get rid of already visited parent in the adj_nodes
  if(parent_node != -1) {
    adj_nodes <- adj_nodes[-which(adj_nodes == parent_node)];
  }
  # Do nothing more if current node is a leaf of a visited branch
  if(length(adj_nodes) > 0) {
    for(n in 1:length(adj_nodes)) {
      create_and_fill_arborescence(
        scores_tab, scores_mat, gene_exp_mat,
        node = adj_nodes[n],
        parent_node = node, parent_struct = node_struct, parent_dir = parent_dir,
        tree_info, tree, sampling_ix);
    }
  }
}

create_and_fill_trajectories <- function(scores_tab, scores_mat, gene_exp_mat,
                                         tree_info, tree, tree_traj,
                                         sampling_ix, dir) {
  for(i in 1:length(tree_traj$paths)) {
    nodes <- tree_traj$paths[[i]];
    traj_dir <- create_subdir(dir, paste("trajectory_", i, sep = ""));
    for(j in 1:length(nodes)) {
      dirname <- paste("Pos", j, "_Node", nodes[j], sep = "");
      node_dir <- create_subdir(traj_dir, dirname);

      node_ix <- which(tree_info$pt_nodes == nodes[j]);

      if(nodes[j] %in% names(sampling_ix)) {
        # scatterplot of sampled individuals within the node
        sampled_node_ix <- sampling_ix[[which(names(sampling_ix) == nodes[j])]][, 1];
        if(any(sampled_node_ix)) {
          scatterplot(node_dir, node_ix[sampled_node_ix], scores_tab, gene_exp_mat);
        }
      }
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
#' With this regard, this function
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
#' @param target_ratio number; fraction of associations to be sampled in each branch
#' @param output_directory where to create the arborescence of directories for each cluster
#'
#' @importFrom ElPiGraph.R computeElasticPrincipalTree PlotPG ConstructGraph GetSubGraph PartitionData
#' @importFrom igraph adjacent_vertices vcount V
#' @importFrom stats dist
#' @importFrom data.table data.table fwrite
explore_clusters <- function(scores_tab, partition, gene_exp_mat,
                             sampling_method = "downsampling", target_ratio = .8,
                             output_directory = "/home/charles/test") {

  # Create a matrix of scores (SPEARMAN rho is ignored) from the table of scores
  # -----------------------------------------------------------------
  scores_mat <- data.matrix(scores_tab[, c("MIC (strength)", "MIC-p^2 (nonlinearity)", "MAS (non-monotonicity)",
                                  "MCN (complexity)", "PEARSON p (linear correlation)")]);

  # Identify clusters from partition
  # -----------------------------------------------------------------
  clusters <- sort(unique(partition));

  # Build a robust tree (resulting averaged tree from bootstrapping)
  # that captures the heterogeneity of each cluster
  # -----------------------------------------------------------------
  trees <- construct_tree(scores_mat, partition);

  # Get relationships between branches, nodes and points of each tree
  # -----------------------------------------------------------------
  trees_info <- get_tree_info(scores_mat, partition, trees);

  # Get the root point of each tree
  # -----------------------------------------------------------------
  trees_starting_pt <- get_clust_starting_pt(scores_mat, partition);
  trees_starting_nodes <- trees_starting_pt;
  for(i in 1:length(clusters)) {
    trees_starting_nodes[[i]] <- trees_info[[i]]$pt_nodes[trees_starting_pt[[i]]];
  }

  # Get the trajectories starting from the root of all trees
  trees_traj_info <-
    get_tree_trajectories(scores_mat, partition,
                          trees, trees_starting_nodes, trees_info);

  # Sample a few points within each node of each tree
  # -----------------------------------------------------------------
  trees_sampling <- list();
  for(i in 1:length(clusters)) {
    trees_sampling[[i]] <-
      sample(scores_tab[partition == clusters[i], ],
             partition = trees_info[[i]]$pt_nodes,
             methods = sampling_method, target_ratio = target_ratio);
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

    # 4. Create the arboresence of directories, meanwhile plotting
    # sampled associations in each node directory (core function)
    # ---------------------------------------------------------------
    create_and_fill_arborescence (
      scores_tab[partition == clusters[i], ],
      scores_mat[partition == clusters[i], ],
      gene_exp_mat,
      node = curr_root_node,
      parent_node = -1, parent_struct = -1,
      parent_dir = curr_cluster_dir,
      tree_info = trees_info[[i]],
      tree = trees[[i]],
      sampling_ix = trees_sampling[[i]]);
  }

  # Create the template csv file to help classifying the associations
  # -----------------------------------------------------------------
  create_classification_template(output_directory, trees_info);

  return(trees_info);
}

output_directory = "~/Documents/vignettes_out/";
partition = partition;
gene_exp_mat = gene_exp_mat;
sampling_method = "downsampling";
target_ratio = .8;
scores_tab = signif_scores;

