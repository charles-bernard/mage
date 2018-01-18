# mage R package

#' @title phenograph_clustering
#'
#' @description uses the phenograph algorithm to perform an automated clustering
#' of the table of scores
#'
#' @param x data.table; the table of scores
#' @param ngb number of nearest neighbours
#'
#' @return
#' a vector corresponding to the partition of the table of scores
#'
#' @importFrom Rphenograph Rphenograph
#' @importFrom igraph membership
phenograph_clustering <- function(x, ngb = 50) {
  x <- x[, 3:ncol(x)];
  Rphenograph_out <- Rphenograph(x[, 3:ncol(x)], k = ngb);
  partition <- membership(Rphenograph_out[[2]]);
  names(partition) <- NULL
  return(as.numeric(partition));
}

#' @title get_optimal_k
#'
#' @description a wrapper to \code{fpc::prediction.strength} function. This function is based on
#' consensus clustering using clara algorithm for different k (2 to 10 precisely).
#'
#' The prediction strength is defined according to Tibshirani and Walther (2005) who recommand to
#' choose as optimal number of clusters the largest number of clusters that leads to
#' a prediction strength above a \code{cutoff} of 0.8 or 0.9.
#'
#' @param x data.table; the table of scores
#' @param ix boolean vector of indexes to which x must be reduced (sampling output)
#' @param cutoff cutoff for the prediction strength (recommanded to range from 0.8 up to 0.9)
#' @param nrep nb of times the clustering is performed for each value of k (bootstrapping)
#'
#' @return the optimal number of clusters
#' @importFrom fpc prediction.strength
get_optimal_k <- function(x, ix = NULL, cutoff = 0.8, nrep = 20) {
  x <- x[, 3:ncol(x)];
  if(!is.null(ix)) {
    x <- x[ix, ];
  }

  k <- k <- prediction.strength(
    x, Gmin = 2, Gmax = 10, M = 20, cutoff = cutoff,
    clustermethod = claraCBI, classification = "centroid")$optimalk;
  return(k);
}

#' @title get_relevant_k
#'
#' @description calls NbClust package to return a list of relevant k
#' according to different methodologies
#'
#' @param x data.frame, the table of scores
#' @param kmax maximun number of clusters allowed
#' @param ix boolean vector of indexes to which x must be reduced (sampling output)
#'
#' @return
#' a list of 2 elements
#' \describe{
#'   \item{NbClust_output}{A vector of integer storing the best k returned by each method}
#'   \item{majority_rule_k}{k most often returned}
#' }
#'
#' @importFrom NbClust NbClust
get_relevant_k <- function(x, ix = NULL, kmax = 12) {
  x <- x[, 3:ncol(x)];
  if(!is.null(ix)) {
    x <- x[ix, ];
  }

  cat("  Computing relevant nb k of clusters ...\n");
  relevant_k <- list();

  # Call NbClust
  # ----------------------------------------------------------------------
  pdf(NULL) # avoid NbClust graphical output
  nbclust_out <- NbClust(x,
                         min.nc = 2, max.nc = kmax,
                         distance = "euclidean", method = "kmeans",
                         index = "all");
  junk <- dev.off();

  # Get the k returned by each indice
  # ----------------------------------------------------------------------
  good_ix <- which(nbclust_out$Best.nc[1, ] >= 1 &
                     nbclust_out$Best.nc[1, ] <= kmax);
  best_k <- nbclust_out$Best.nc[1, good_ix];
  best_k_names <- names(best_k);
  relevant_k$`NbClust_output` <- as.integer(best_k);
  names(relevant_k$`NbClust_output` ) <- best_k_names;

  # Majority Rule k
  # ----------------------------------------------------------------------
  t <- table(relevant_k$`NbClust_output`);
  relevant_k$`majority_rule_k` <-  as.integer(names(t)[max(which(t == t[which.max(t)]))]);
  return(relevant_k);
}

#' @title clara_clustering
#'
#' @description uses clara algorithm to partition the table of scores into k clusters
#' around medoids (a more robust version of K-means)
#'
#' @param x data.frame; the table of scores
#' @param k desired number of clusters
#'
#' @return
#' a vector corresponding to the partition of the table of scores
#'
#' @importFrom fpc claraCBI
clara_clustering <- function(x, k) {
  # Find the best partition for k clusters using clara algorithm
  # ----------------------------------------------------------------------
  x <- x[, 3:ncol(x)];
  res.clara <- claraCBI(x, k, usepam = TRUE);
  return(res.clara$partition);
}



