# mage R package

# This file contains a big function to sample a table of scores or clusters of it
# according to different methods ((down|low|up) sampling, (centroids|medoids) selection)

# Copied from spade R package (Multidimensional kernel density estimation)
#' @useDynLib mage compute_density_
compute_density <- function(x, kernel_mult = 5.0, apprx_mult = 1.5, med_samples = 2000) {
  .Call("compute_density_", t(x), kernel_mult, apprx_mult, med_samples);
}

get_boundary <- function(density, target_nb) {

  density_s <- sort(density, decreasing = TRUE);
  cdf <- rev(cumsum(1.0 / density_s));
  boundary <- target_nb / cdf[1];

  # Default solution if target density smaller than any present
  if(boundary > density_s[length(density_s)]) {
    targets <- (target_nb - 1:length(density_s)) / cdf;
    return(targets[which.min(targets - rev(density_s) > 0)]);
  }

  return(boundary);
}


#' @title sample
#'
#' @description
#' Sample the table of scores or clusters of it
#' according to different methods of selection
#'
#' @param x data.frame, a table of scores
#' @param partition vector containing the partition of table \code{x}.
#'
#' If \code{partition} is passed to this function, each cluster will be sampled instead of the whole table.
#'
#' Each element of the partition vector must be a cluster ID and length(partition) must be equal to nrow(x).
#'
#'
#' @param methods method or vector of methods to be used for the sampling.
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
#' See Details to have the description of the methods
#'
#' @param target_ratio percentage of gene pairs to be retained in \code{x}
#' (if partition is NULL) or in each cluster of \code{x} (if partition is NOT NULL)
#'
#' argument \code{target_ratio} will be ignored if argument \code{target_nb} is provided
#'
#' @param target_nb number of gene pairs to be retained in \code{x} (if partition is NULL)
#' or in each cluster of \code{x} (if partition is NOT NULL)
#'
#' @importFrom stats dist median quantile runif
sample <- function(x,
                   partition = NULL,
                   methods = "downsampling",
                   target_ratio = .1,
                   target_nb = NULL) {

  x <- data.matrix(x[, 3:ncol(x)]);
  sampling_out <- list();

  if(is.null(methods)) {
    stop ("Methods argument is NULL");
  }

  # Check whether method name is correct
  # --------------------------------------------------------------------
  for(curr_method in methods){
    if(!curr_method %in% c("downsampling", "lowsampling", "upsampling",
                           "centroids", "medoids")) {
      stop("at least one method name is not expected");
    }
  }

  # Define target_nb if not already defined
  # --------------------------------------------------------------------
  if(is.null(target_nb)) {
    if(is.null(target_ratio)) {
      stop("A target ratio or nb must be provided");
    } else {
      computed_target_nb <- round(nrow(x) * target_ratio);
    }
  } else {
    computed_target_nb <- target_nb;
    if(target_nb > nrow(x)) {
      computed_target_nb <- nrow(x);
    }
  }

  # If no partition is provided ->
  #   the function follows the instructions
  # ----------------------------------------------------------------------
  if(is.null(partition)) {

    # Initialize table of boolean indexes to be returned
    # --------------------------------------------------------------------
    sampling_table <- data.frame(matrix(FALSE, nrow = nrow(x), ncol = length(methods)));
    colnames(sampling_table) <- methods;

    # Compute density if there are density dependent sampling methods
    # --------------------------------------------------------------------
    for(curr_method in methods) {
      if(curr_method %in% c("downsampling", "lowsampling", "upsampling")) {
          density <- compute_density(x);
          random_proba <- runif(length(density));
          break;
      }
    }

    for(curr_method in methods){
      curr_col <- which(colnames(sampling_table) == curr_method);

      # Density Dependent Sampling methods
      # ------------------------------------------------------------------
      if(curr_method %in% c("downsampling", "lowsampling", "upsampling")) {
        curr_density <- density;
        curr_proba <- random_proba;
        considered_ix <- 1:length(density);

        if(curr_method == "upsampling") {
          # Invert vector of densities to maximize high density individuals
          curr_density <- max(density) - density + min(density);

        } else if(curr_method == "lowsampling") {
          # Rarest individuals need to be identified because they will be eventually filtered out
          rarest_density_thresh <- quantile(density, .01) + ceiling(length(density)*.001);
          rarest_ix <- which(density < rarest_density_thresh);
          considered_ix <- considered_ix[-rarest_ix];
          curr_proba[rarest_ix] <- max(density) + 1; # This makes sure that rarest_idx wont be retained
        }

        # Sampling
        curr_boundary <- get_boundary(curr_density[considered_ix], computed_target_nb);
        sampling_table[, curr_col] <- curr_proba < (curr_boundary / curr_density);
      }

      # Centroids and Medois Selection
      # ------------------------------------------------------------------
      if(curr_method %in% c("centroids", "medoids")) {
        if(curr_method == "centroids") {
          reference_coordinates <- colMeans(x);
        } else {
          reference_coordinates <- apply(x, 2, median);
        }

        # Get nearest neighbor from (centr|med)oïd
        alldist <- as.matrix(dist(rbind(reference_coordinates, x)));
        mindist_relative_to_ref <- sort(alldist[1,-1], index.return = T)$ix;
        sampling_table[mindist_relative_to_ref[1:computed_target_nb], curr_col] = TRUE;
      }
    }
    sampling_out$`Whole Table` <- sampling_table;
    return(sampling_out);
  }

  # Else if a partition is provided ->
  #   Recursivity with x = current_custer_table and partition = NULL
  # ----------------------------------------------------------------------
  if(!is.null(partition)) {
    clusters <- unique(partition);

    for(i in 1:length(clusters)) {
      cluster <- clusters[i];
      cluster_ix <- which(partition == cluster);
      sampling_out[[i]] <- sample(x[cluster_ix, ], partition = NULL,
                                  methods = methods,
                                  target_ratio = target_ratio, target_nb = target_nb)[[1]];
      names(sampling_out)[i] <- cluster;
    }

    return(sampling_out);
  }
}

