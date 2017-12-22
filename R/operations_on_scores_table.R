# mage R package

# This file regroups functions to:
#   1. Plots the scatterplot of the scores (to visualize their relationships)
#   2. Filter the matrix based on a threshold of significance
#   3. Standardize the matrix (using zscore)


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


#' @title filter_scores
#'
#' @description
#' Filters a table of scores based on a threshold of significance defined either by a MIC score, or by a p-value
#'
#' @param x data.frame: the table of scores
#' @param on variable used for the threshold (accepts only \code{"MIC"} or \code{"pval"})
#' @param pval numeric vector of p-values corresponding to the MIC scores in \code{x} (only if \code{on = "pval"})
#' @param thresh value of the threshold
#'
#' @return
#' Returns a filtered table of scores
#'
#' @examples
#' \dontrun{
#' scores_tab <- compute_scores(my_gene_exp_matrix, n.cores = 6)
#' pvalues <- assign_pval(scores_tab$`MIC (strength)`, nb_cells = 96)
#' signif_scores_tab <- filter_scores(scores_tab, on = "pval", 
#'                                    pval = pvalues, threshold = 0.05)
#' }                                    
#'
#' @details
#' if \code{on = "MIC"}, \code{filter_scores} will retain any association whose MIC >= \code{threshold}
#' 
#' else if \code{on = "pval"}, \code{filter_scores} will retain any association whose p-value <= \code{threshold}
#'   
filter_scores <- function(x,
                   on = 'MIC',
                   pval = NULL,
                   thresh = 0.4) {
  # To do: check all conditions
  
  if(on == "MIC") {
    if(thresh < 0 || thresh > 1) {
      stop("A MIC score is defined in [0;1]");
    }
    retained_ix <- which(x[, 3] >= thresh);
  } else if(on == "pval") {
    if(is.null(pval)) {
      stop("argument pval missing");
    }
    if(length(pval) != nrow(x)) {
      stop("length(pval) must be equal to nrow(x)");
    }
    retained_ix <- which(pval <= thresh);
  } else {
    stop("on argument must be set to either \"MIC\" or \"pval\"");
  }
  
  return(x[retained_ix, ]);
}

#' @title standardize_scores
#'
#' @description
#' Standardize each column of the table of scores via the Zscore: \deqn{x = ( x - mean(column) ) / ( std(column) )}
#'
#' @param x data.frame, table of scores
#' 
#' @return
#' returns a standardized table of scores
standardize_scores <- function(x) {
  std_x <- data.frame(scale(x[, 3:ncol(x)]));
  
  possible_NA_col <- apply(std_x, 2, function(i) any(is.na(i)));
  if(length(which(possible_NA_col == TRUE)) > 0) {
    std_x <- std_x[, -(which(possible_NA_col == TRUE))];
  }
  
  return(data.frame(x[, 1:2], std_x));
}


# file <- "/media/charles/Seagate Expansion Drive/Curie/Analyses/Analyses_CDD/MAGE/A471/794_remaining_mvg/794vg_Mage_out/1-Association_Scores/0-scores.csv";
# data <- fread(file, sep = ",", header = T);
# MIC_scores <- data[,3][[1]];
# nb_cells <- 96;


