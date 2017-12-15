# mage R package

# This file regroups functions to:
#   1. Compute the table of scores (both asso. and charac. scores)
#   2. Filter the matrix based on a MIC threshold of significance
#   3. Standardize the matrix (using zscore)

#' @title compute_scores
#'
#' @description
#' Computes the table of association and characterization scores
#'
#' @details Structure of the returned table of scores
#' \itemize{
#'   \item{column 1:2} {- the two genes whose assocation is being tested}
#'   \item{column 3} {- MIC score, which is the strength of the association between the 2 genes}
#'   \item{column 4:n} {- characterization scores, informing about the "shape" of the association}
#' }
#'
#' @param x a gene expression matrix (row.names corresponding to cells, col.names to genes)
#' @param n.cores nb of cores to be used for parallel processing
#'
#'
compute_scores <- function(x,
                           n.cores = 1) {

  n_genes <- nrow(x);

  # Get unique pairs
  pairs <- t(combn(nrow(x), 2));
  n_pairs <- nrow(pairs);
  pairs_name <- paste(rownames(x)[pairs[, 1]],
                      rownames(x)[pairs[, 2]],
                      sep = "-");
  pairs_ix <- 1:n_pairs;

  # Compute MINE scores
  cat(" * Compute the MINE scores for all pairs of variables ...\n");
  cat("     this may take quite some time ...\n");
  MINE_sym_mat <- mine(t(x), master = 1:n_genes,
                     n.cores = n.cores,
                     alpha = 0.6, C = 15);
  n_scores <- length(MINE_sym_mat);
  scores_name <- names(MINE_sym_mat);

  # Turn symetric MINE_matrix to a matrix
  # whose each row is a unique association
  MINE_mat <- matrix(nrow = n_pairs, ncol = n_scores);
  colnames(MINE_mat) <- scores_name;
  for(i in 1:n_scores) {
    MINE_mat[, i] <- unlist(vapply(pairs_ix,
                            function(ix)
                              as.list(MINE_sym_mat[[i]][pairs[ix, 1], pairs[ix, 2]]),
                            FUN.VALUE = c(numeric)));
  }

  # Compute Pearson and Spearman Correlations as well
  cat(" * Compute the Pearson Correlation for all pairs of variables ...\n");
  PEARSON <- vapply(pairs_ix,
                    function(ix) as.list(cor(x[pairs[ix, 1], ],
                                             x[pairs[ix, 2], ],
                                             method = "pearson")),
                    FUN.VALUE = c(numeric));
  cat(" * Compute the Spearman Correlation for all pairs of variables ...\n");
  SPEARMAN <- vapply(pairs_ix,
                     function(ix) as.list(cor(x[pairs[ix, 1], ],
                                              x[pairs[ix, 2], ],
                                              method = "spearman")),
                     FUN.VALUE = c(numeric));

  # Sort MINE_mat by decreasing order of MIC Score
  MINE_table <- data.table(MINE_mat);
  sort_ix <- sort(unlist(MINE_table[,1]),
                  decreasing = T, index.return = T)$ix;

  # Write final table, assign customed column names
  final_table <- data.table(`X gene` = rownames(x)[pairs[sort_ix, 1]],
                            `Y gene` = rownames(x)[pairs[sort_ix, 2]],
                            `MIC (strength)` = MINE_table[sort_ix]$`MIC`,
                            `MIC-p^2 (nonlinearity)` = MINE_table[sort_ix]$`MICR2`,
                            `MAS (non-monotonicity)` = MINE_table[sort_ix]$`MAS`,
                            `MEV (functionality)` = MINE_table[sort_ix]$`MEV`,
                            `MCN (complexity)` = MINE_table[sort_ix]$`MCN`,
                            `PEARSON p (linear correlation)` = PEARSON_results[sort_ix],
                            `SPEARMAN rho (rank correlation)` = SPEARMAN_results[sort_ix]);

  return(final_table);
}

# file = "/home/cbernard/Documents/Datasets/scRnaSeq/Datasets/CAF_S1/Biopsy/A471/most_variant_genes/1000mvg/794_variant_genes_data.csv";
# x = data.matrix(read.table(file, sep = ",", header = T, row.names = 1));
