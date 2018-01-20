#' Single cell gene expression matrix of Peripheral Blood Mononuclear Cells (PBMC)
#'
#' A dataset containing the raw count for 230 genes in 80 cells.
#'
#' This dataset is a reduced version of the PBMC dataset (2700 cells)
#' freely available from 10X genomics at:
#'
#' \url{https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz}
#'
#' This dataset has been shorten by the Seurat R package to serve
#' as canonical input for the code examples in the documentation of the functions
#'
#' The Seurat R object corresponding to this example dataset can be found
#' at:
#'
#' \url{https://github.com/satijalab/seurat/blob/master/data/pbmc_small.rda}
#'
#' This matrix has been produced by loading the Seurat R object and by
#' extracting the raw count sparse matrix as follow:
#'
#' \code{pbmc_samall_raw_data <- data.matrix(pbmc_small@@raw.data);}
#'
#' @format A matrix with 230 rows (genes) and 80 columns (cells)
"pbmc_small_raw_data"
