# mage R package

# This file regroups functions to:
#   1. Plots the scatterplot of the scores (to visualize their relationships)


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
#' @param pdf_filename path to the pdf file (optional)
#'
plot_scores_relationships <- function(x, pdf_filename = NULL) {
  scores_mat <- data.matrix(x[, 3:ncol(x)]);
  n_scores <- ncol(scores_mat);

  if(!is.null(pdf_filename)) {
    grDevices::pdf(pdf_filename, width = 10, height = 10);
  }

  graphics::par(mfrow = c(n_scores, n_scores),
                mai = c(.2, .2, .2, .2),
                oma = c(2, 2, .2, .2));

  for(i in 1:n_scores) {
    for(j in 1:n_scores) {
      if(j <= i) {
        print(sprintf("Plotting %s vs %s ...", colnames(scores_mat)[i], colnames(scores_mat)[j]));
        graphics::par(mfg = c(i, j));
        graphics::plot(scores_mat[, j], scores_mat[, i], col = "blue", pch = 16, cex = .5, ann = F);
        if(j == 1)
          graphics::mtext(colnames(scores_mat)[i], side = 2, cex = 0.7, padj = -4);
        if(i == j)
          graphics::mtext(colnames(scores_mat)[j], side = 3, cex = 0.7, padj = -1);
      }
    }
  }

  if(!is.null(pdf_filename)) {
    junk <- grDevices::dev.off();
  }
}


