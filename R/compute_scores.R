# mage R package

# This file contains a single function to:
#   1. Compute the table of scores (both asso. and charac. scores)

#' @title compute_scores
#'
#' @description
#' Computes the table of association and characterization scores
#'
#' @param x matrix: a gene expression matrix (rownames corresponding to cells and colnames to genes)
#' @param n.cores integer: nb of cores to be used for parallel processing
#'
#' @return
#' Returns a table with the following columns:
#' \describe{
#'   \item{1:2}{Names of the two genes whose assocation is being assessed}
#'   \item{3}{MIC score, \emph{i.e} the strength of the association between the 2 genes}
#'   \item{4:n}{Characterization scores, informing about the "shape" of the association}
#' }
#'
compute_scores <- function(x,
                           n.cores = 1) {

  n_genes <- nrow(x);

  # Get unique pairs
  # ----------------------------------------------------------------------
  pairs <- t(utils::combn(nrow(x), 2));
  n_pairs <- nrow(pairs);
  pairs_name <- paste(rownames(x)[pairs[, 1]],
                      rownames(x)[pairs[, 2]],
                      sep = "-");
  pairs_ix <- 1:n_pairs;

  # Compute MINE scores
  # ----------------------------------------------------------------------
  cat(" * Compute the MINE scores for all pairs of variables ...\n");
  cat("     this may take quite some time ...\n");
  MINE_sym_mats <- minerva::mine(t(x), master = 1:n_genes,
                                 n.cores = n.cores,
                                 alpha = 0.6, C = 15);
  n_scores <- length(MINE_sym_mats);
  scores_name <- names(MINE_sym_mats);

  # Turn the symetric matrices into a global table
  # whose each row is a unique association
  MINE_mat <- matrix(nrow = n_pairs, ncol = n_scores);
  colnames(MINE_mat) <- scores_name;
  for(i in 1:n_scores) {
    MINE_mat[, i] <- unlist(vapply(pairs_ix,
                            function(ix)
                              as.list(MINE_sym_mats[[i]][pairs[ix, 1], pairs[ix, 2]]),
                            FUN.VALUE = c(numeric)));
  }

  # Compute Pearson and Spearman Correlations as well
  # ----------------------------------------------------------------------
  cat(" * Compute the Pearson Correlation for all pairs of variables ...\n");
  PEARSON_cor <- vapply(pairs_ix,
                        function(ix) as.list(stats::cor(x[pairs[ix, 1], ],
                                                    x[pairs[ix, 2], ],
                                                    method = "pearson")),
                        FUN.VALUE = c(numeric));
  cat(" * Compute the Spearman Correlation for all pairs of variables ...\n");
  SPEARMAN_cor <- vapply(pairs_ix,
                         function(ix) as.list(stats::cor(x[pairs[ix, 1], ],
                                                     x[pairs[ix, 2], ],
                                                     method = "spearman")),
                         FUN.VALUE = c(numeric));

  # Sort MINE_mat by decreasing order of MIC Score
  # ----------------------------------------------------------------------
  MINE_table <- data.frame(MINE_mat);
  sort_ix <- sort(unlist(MINE_table[,1]),
                  decreasing = T, index.return = T)$ix;

  # Write final table, assign customed column names
  # ----------------------------------------------------------------------
  final_table <- data.frame(`X gene` = rownames(x)[pairs[sort_ix, 1]],
                            `Y gene` = rownames(x)[pairs[sort_ix, 2]],
                            `MIC (strength)` = MINE_table[sort_ix]$`MIC`,
                            `MIC-p^2 (nonlinearity)` = MINE_table[sort_ix]$`MICR2`,
                            `MAS (non-monotonicity)` = MINE_table[sort_ix]$`MAS`,
                            `MEV (functionality)` = MINE_table[sort_ix]$`MEV`,
                            `MCN (complexity)` = MINE_table[sort_ix]$`MCN`,
                            `PEARSON p (linear correlation)` = PEARSON_cor[sort_ix],
                            `SPEARMAN rho (rank correlation)` = SPEARMAN_cor[sort_ix]);

  return(final_table);
}
