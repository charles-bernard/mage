# mage R package

# This file contains a big function to sample a table of scores or clusters of it
# according to different methods ((down|low|up) sampling, (centroids|medoids) selection)

# Copy from spade R package
#' @useDynLib mage compute_density_
compute_density <- function(x, kernel_mult=5.0, apprx_mult=1.5, med_samples=2000)
  .Call("compute_density_",t(x),kernel_mult,apprx_mult,med_samples)

# #' @title sample
# #'
# #' @description
# #' Computes the table of association and characterization scores
# #'
# sample <- function(x,
#                    partition = NULL,
#                    methods = "downsampling",
#                    target_nb = NULL,
#                    target_ratio = .1) {
# }



