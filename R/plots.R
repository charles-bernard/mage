# mage R package

# This file regroups functions to:
#   1. Plots the scatterplot of the scores (to visualize their relationships)
#   2. Show heatmap of scores for a provided partition
#   3. Show individuals in a space reduced to PCs

#' @title plot_scores_relationships
#'
#' @description
#' Produce the scatterplots of each pair of scores and display them
#' as a symetric matrix of plots.
#'
#' This function is meant to visualize
#' the kind of relationships that exist between the different scores.
#'
#' @param x data.frame: the table of scores
#' @param pdf_filename path to a pdf file which must store the plot (optional)
#'
#' @importFrom graphics par plot mtext
#' @importFrom grDevices pdf dev.off
plot_scores_relationships <- function(x, pdf_filename = NULL) {
  x <- data.matrix(x[, 3:ncol(x)]);
  n_scores <- ncol(x);

  if(!is.null(pdf_filename)) {
    pdf(pdf_filename, width = 10, height = 10);
  }

  par(mfrow = c(n_scores, n_scores),
      mai = c(.2, .2, .2, .2),
      oma = c(2, 2, .2, .2));

  for(i in 1:n_scores) {
    for(j in 1:n_scores) {
      if(j <= i) {
        print(sprintf("Plotting %s vs %s ...", colnames(x)[i], colnames(x)[j]));
        par(mfg = c(i, j));
        plot(x[, j], x[, i], col = "blue", pch = 16, cex = .5, ann = F);
        if(j == 1)
          mtext(colnames(x)[i], side = 2, cex = 0.7, padj = -4);
        if(i == j)
          mtext(colnames(x)[j], side = 3, cex = 0.7, padj = -1);
      }
    }
  }

  if(!is.null(pdf_filename)) {
    junk <- dev.off();
  }
}

#' @title show_heatmap
#'
#' @description
#' Show the heatmap of scores for each cluster in the partition
#'
#' @param x data.frame; the standardize table of scores (consider using \code{standardize_scores} function)
#' @param partition vector containing the partition of the table \code{x}.
#' @param ix boolean vector of indexes to which the heatmap must be reduced (\code{sample} output).
#' A hundred representative individuals by cluster (downsampling) should be enough to get a glimpse
#' of what the cluster looks like in term of profile of scores ...
#'
#' \code{ix} is a recommanded argument for a large table of score (greater than 1000 individuals)
#' @param pdf_filename path to a pdf file which must store the plot (optional)
#'
#' @import methods
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics par plot mtext
#' @importFrom grDevices pdf dev.off
#' @import ComplexHeatmap
show_heatmap <- function(x, partition, ix = NULL, pdf_filename = NULL) {
  # Some basic controls about the provided partition
  if(is.null(partition)) {
    stop("a partition must be provided")
  } else {
    if(length(partition) != nrow(x)) {
      stop("length of the partition must be equal to nrow(x)");
    }
  }

  x <- x[, 3:ncol(x)];
  if(!is.null(ix)) {
    x <- x[ix, ];
    partition <- partition[ix];
  }

  colorscale <- c(brewer.pal(n = 9, name = "Set1"),
                  brewer.pal(n = 8, name = "Set2"),
                  brewer.pal(n = 8, name = "Accent"));

  cat(" * Plotting the heatmap ...\n");

  # Heatmap plot
  # ----------------------------------------------------------------------
  if(!is.null(pdf_filename)) {
    pdf(pdf_filename);
  }
  par(mar = c(3, 4, 2, 4) + 0.1,
      oma = c(3, 3, 1, 3))
  if(requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    hm <- Heatmap(x, cluster_columns = FALSE,
                  clustering_distance_rows = "pearson",
                  clustering_method_rows = "ward.D2",
                  split = partition,
                  name = "Standardized Scores");
    hm <- hm +
      Heatmap(partition,
              name = "Partition", width = grid::unit(15, "mm"),
              col = colorscale[1:length(unique(partition))],
              cluster_columns = FALSE, cluster_rows = FALSE,
              heatmap_legend_param = list(color_bar = "discrete"));
    print(hm);
  } else {
    stop("package 'ComplexHeatmap' not found, please install this package")
  }
  if(!is.null(pdf_filename)) {
    junk <- dev.off();
  }
}

# Internal function: Plot report of the PCA
# ----------------------------------------------------------------------
pca_report <- function(pca_out) {
  scores <- get_pca_var(pca_out);
  p <- list()

  # PLOT POURCENTAGE OF EXPLAINED VARIANCE PER PCA DIMENSION
  p$percentage_of_explained_variance_per_PC <-
    fviz_eig(pca_out,
             addlabels = TRUE,
             main = "Percentage of Explained Variance per PCA Dimension",
             font.main = 20,
             font.x = 15,
             font.y = 15,
             font.tickslab = 15);

  # a little trick to store the corrplot inside a variable
  store_corrplot <- function(x) {
    cp <- function() {
      par(oma=c(5,5,5,5));
      corrplot(x,
               is.corr=FALSE,
               title = "Quality of score representation (cos2) per PCA Dimension",
               mar = c(5,4,6,2) + 0.1);
    }
    list(show = cp);
  }

  # QUALITY OF VARIABLE REPRESENTATION PER PCA DIM
  p$quality_of_variable_representation_per_PC$message <-
    "Add suffix PCA_information$quality_of_variable_representation_per_PC$show() to visualize the corrplot";
  p$quality_of_variable_representation_per_PC$show <- store_corrplot(scores$cos2)$show;

  # CONTRIBUTION OF VARIABLES TO EACH PCA DIMENSION
  p$variables_contribution_to_PCs$message <-
    "Add_suffix PCA_information$variables_contribution_to_PCs$show() to visualize the corrplot";
  p$variables_contribution_to_PCs$show <- store_corrplot(scores$contrib)$show;

  class(p$quality_of_variable_representation_per_PC) <-
    class(p$variables_contribution_to_PCs) <- "corrplot";
  return(p);
}

# Internal function: turn boolean indexes into vector of sampling classes
# ----------------------------------------------------------------------
ix_to_factors <- function(ix, partition) {
  n <- length(ix);
  length_factors <- 0;
  # nb clusters * nrow by cluster gives the length of the vector
  for(i in 1:n) {
    length_factors <- length_factors + nrow(ix[[i]]);
  }
  if(is.null(partition)) {
    partition <- rep(names(ix)[1], length_factors);
  }

  factors <- rep("filtered", length_factors);
  for(i in 1:n) {
    # Get indexes of the current cluster
    curr_ix <- which(partition == names(ix)[i]);
    # Assign "method class" to retained individuals within clusters
    for(j in 1:ncol(ix[[i]])) {
      factors[curr_ix][ix[[i]][, j]] <- colnames(ix[[i]])[j];
    }
  }

  return(factors);
}

# Internal function: Plot individuals on PCs space
# ----------------------------------------------------------------------
plot_individuals_on_pc <- function(pca_out, my_habillage, PCs = 1:5, ell = FALSE) {
  PCs <- sort(PCs);
  colorscale <- c("grey",
                  brewer.pal(n = 9, name = "Set1"),
                  brewer.pal(n = 8, name = "Set2"),
                  brewer.pal(n = 8, name = "Accent"));

  p <- list(); k  <- 1;
  for(i in PCs) {
    for(j in PCs) {
      if(i < j) {
        if(i == 1 && j == 2) {
          ell = ell;
        } else {
          ell = FALSE;
        }
        if(is.null(my_habillage)) {
          my_habillage <- "none";
          my_palette <- NULL;
        } else {
          my_palette <- colorscale[1:length(unique(my_habillage))];
          my_habillage <- as.factor(my_habillage);
        }
        p[[k]] <-
          fviz_pca_biplot(pca_out,
                          axes = c(i,j),
                          repel = F,
                          font.main = 20,
                          font.x = 15,
                          font.y = 15,
                          font.tickslab = 15,
                          palette = my_palette,
                          # individuals
                          point.size = 0.5,
                          habillage = my_habillage,
                          alpha.ind = 0.7,
                          # variables
                          label = c("var"),
                          col.var = "black",
                          addEllipses = ell);
        names(p)[k] <- paste("PC", i, "vs_PC", j, sep = "_");
        k <- k + 1;
      }
    }
  }

  return(p);
}

#' @title show_reduced_space
#' @description
#' Reduce the space of the scores to its first principal components
#' in order to display on it the associated gene pairs.
#'
#' @param x data.frame; the table of scores
#' @param x_ixs (optional) list containing a unique n-by-m matrix of boolean indexes
#' resulting from the \code{sample} function
#' where \code{n} is the number of individuals in \code{x}
#' and \code{m} is the number of methods used for the sampling
#' @param partition  (optional) vector; partition of \code{x} into k clusters
#' @param partition_ixs (optional) list of n-by-m matrices of boolean indexes
#' resulting from \code{sample} function
#' where \code{n} is the number of individuals in a given cluster of the \code{partition}
#' and \code{m} is is the number of methods used for the sampling
#'
#' @param PCs numeric vector of Principal Compenents to be considered (max = 5)
#' @param output_dir (optional) path to a directory where to store all the plots
#'
#' @details
#' the function \code{show_reduced_space} returns a graphical object
#' composed of 2 up to 5 elements:
#' \describe{
#'   \item{PCA_information}{The graphical informations about the PCA}
#'   \item{population}{Show all the gene pairs in the reduced space}
#'   \item{population_samples}{Show the sampled gene pairs within the population
#'   (\code{ixs} must be provided)}
#'   \item{clusters}{Show the different clusters of the population
#'   (\code{partition} must be provided)}
#'   \item{clusters_samples}{Show the different sampled gene pairs within each clusters
#'   (\code{partition_ixs} must be provided)}
#' }
#'
#' @return
#' Returns a list of plots
#'
#' @importFrom graphics par plot
#' @importFrom grDevices pdf dev.off
#' @importFrom corrplot corrplot
#' @importFrom FactoMineR PCA
#' @importFrom factoextra get_pca_var fviz_eig fviz_pca_var fviz_pca_biplot
show_reduced_space <- function(x, x_ixs = NULL,
                               partition = NULL, partition_ixs = NULL,
                               PCs = 1:4,
                               output_dir = NULL) {
  x <- data.matrix(x[, 3:ncol(x)]);
  pca_out <- PCA(x, graph = FALSE);

  reduced_space <- list();

  # PCA information
  reduced_space$PCA_information <- pca_report(pca_out);

  # Population
  reduced_space$population <-
    plot_individuals_on_pc(pca_out,
                           my_habillage = NULL,
                           PCs);

  # Sampled individuals within population
  if(!is.null(x_ixs)) {
    habillage = ix_to_factors(x_ixs, partition = NULL);
    reduced_space$population_samples <-
      plot_individuals_on_pc(pca_out,
                             my_habillage = habillage,
                             PCs);
  }

  # Clusters
  if(!is.null(partition)) {
    reduced_space$clusters <-
      plot_individuals_on_pc(pca_out,
                             my_habillage = partition,
                             PCs, ell = TRUE)
  }

  # Sampled individuals within clusters
  if(!is.null(partition_ixs)) {
    habillage = ix_to_factors(partition_ixs, partition = partition);
    reduced_space$clusters_samples <-
      plot_individuals_on_pc(pca_out,
                             my_habillage = habillage,
                             PCs);
  }


  cat("NB: each element of the returned list is a plot\n");
  cat("Structure of the returned list:\n\n")
  print(lapply(reduced_space, names));

  # Print the plots in pdf files
  # ----------------------------------------------------------------------
  if(!is.null(output_dir)) {
    cat("Storing the plots in pdf files inside the output directory ...");
    # PCA info
    # ----------------------------------------------------------------------
    pca_info_file <- file.path(output_dir, "PCA_information.pdf");
    pdf(pca_info_file, width = 10, height = 10);
    plot(reduced_space$PCA_information$percentage_of_explained_variance_per_PC);
    reduced_space$PCA_information$quality_of_variable_representation_per_PC$show();
    reduced_space$PCA_information$variables_contribution_to_PCs$show();
    junk <- dev.off();

    # Population
    # ----------------------------------------------------------------------
    population_file <- file.path(output_dir, "population.pdf");
    pdf(population_file, width = 10, height = 10);
    for(i in 1:length(reduced_space$population)) {
      plot(reduced_space$population[[i]]);
    }
    junk <- dev.off();

    if(!is.null(reduced_space$population_samples)) {
      population_samples_file <- file.path(output_dir, "population_samples.pdf");
      pdf(population_samples_file, width = 10, height = 10);
      for(i in 1:length(reduced_space$population_samples)) {
        plot(reduced_space$population_samples[[i]]);
      }
      junk <- dev.off();
    }

    if(!is.null(reduced_space$clusters)) {
      clusters_file <- file.path(output_dir, "clusters.pdf");
      pdf(clusters_file, width = 10, height = 10);
      for(i in 1:length(reduced_space$clusters)) {
        plot(reduced_space$clusters[[i]]);
      }
      junk <- dev.off();
    }

    if(!is.null(reduced_space$clusters_samples)) {
      clusters_samples_file <- file.path(output_dir, "clusters_samples.pdf");
      pdf(clusters_samples_file, width = 10, height = 10);
      for(i in 1:length(reduced_space$clusters_samples)) {
        plot(reduced_space$clusters_samples[[i]]);
      }
      junk <- dev.off();
    }
  }

  return(reduced_space);
}

# internal function; create subdir
# ----------------------------------------------------------------------
create_subdir <- function(main_dir, subdirname) {
  subdir = file.path(main_dir, subdirname);
  if(file.exists(subdir))
    unlink(subdir, recursive = TRUE);
  dir.create(subdir);
  return(subdir);
}

# internal function: scatterplot
# ----------------------------------------------------------------------
scatterplot <- function(output_dir, ix, x, y) {
  for(i in 1:length(ix)) {
    X_gene <- x[ix[i], 1];
    Y_gene <- x[ix[i], 2];
    X_data <- c(y[as.character(X_gene), ]);
    Y_data <- c(y[as.character(Y_gene), ]);
    plot_file <- file.path(output_dir, paste(i, "-", X_gene, "_vs_", Y_gene, ".svg", sep = ""));
    svg(plot_file, width = 20, height = 7);
    par(mfcol = c(1,3));
    plot(as.numeric(X_data),
         as.numeric(Y_data),
         pch = 19,
         xlab = X_gene,
         ylab = Y_gene,
         main = paste(X_gene, "Expression vs", Y_gene, "Expression"),
         sub = paste("MIC =", round(x[ix[i], 3], 3)));
    plot(as.numeric(Y_data),
         as.numeric(X_data),
         pch = 19,
         xlab = Y_gene,
         ylab = X_gene,
         main = paste(Y_gene, "Expression vs", X_gene, "Expression"),
         sub = paste("MIC =", round(x[ix[i], 3], 3)));
    plot(log2(as.numeric(X_data)+1),
         log2(as.numeric(Y_data)+1),
         pch = 19,
         xlab = X_gene,
         ylab = Y_gene,
         main = paste("log2(", X_gene, "Expression ) vs log2(", Y_gene, "Expression )"),
         sub = paste("MIC =", round(x[ix[i], 3], 3)));
    junk <- dev.off();
  }
}

#' @title show_associations
#' @description
#' Produces the scatterplot of each sampled association between 2 genes inside
#' each cluster of a partition
#'
#' @param x data.frame; the table of scores
#' @param y matrix; the gene expression matrix
#' (rownames corresponding to genes and colnames to cells)
#' @param partition vector; partition of \code{x} into k clusters
#' @param partition_ixs list of n-by-m matrices of boolean indexes
#' resulting from \code{sample} function
#' where \code{n} is the number of individuals in a given cluster of the \code{partition}
#' and \code{m} is is the number of methods used for the sampling
#' @param output_dir path to a directory where to store all the plots
#'
#' @importFrom grDevices svg dev.off
#' @importFrom graphics par plot
show_associations <- function(x, y,
                              partition,
                              partition_ixs,
                              output_dir) {
  if(is.null(partition)) {
    stop("A partition must be provided");
  }
  if(is.null(partition_ixs)) {
    stop("Indexes of sampled individuals in each cluster of the partition must be provided");
  }
  if(is.null(y)) {
    stop("A gene expression matrix must be provided");
  }
  if(is.null(output_dir)) {
    stop("An output_dir must be provided");
  }

  n <- length(partition_ixs);
  # for each cluster
  for(i in 1:n) {
    clust <- names(partition_ixs)[i];
    clust_dirname <- paste("cluster", clust, sep = "_");
    clust_dir <- create_subdir(output_dir, clust_dirname);
    clust_ix <- which(partition == clust);
    cat(sprintf("Plotting sampled associations within cluster %s\n", clust));
    # for each sampling method
    for(j in 1:ncol(partition_ixs[[i]])) {
      method_dirname <- colnames(partition_ixs[[i]])[j];
      method_dir <- create_subdir(clust_dir, method_dirname);
      sampled_ix <- clust_ix[partition_ixs[[i]][, j]];
      scatterplot(method_dir, sampled_ix, x, y);
    }
  }
}
