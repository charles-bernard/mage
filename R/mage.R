#' @title
#' mage: Mine Associated Gene Expressions (from single-cell Rna Seq data)
#' 
#' @description
#' The purpose of mage is to capture, and more essentially to classify,
#' any association between 2 gene expressions in a Single Cell Rna Seq
#' Dataset.
#' 
#' @details
#' \describe{
#'     \item{acronym}{Mine Associated Gene Expressions}
#'     \item{relies on}{The Mine statistics to measure the strength of association between 2 gene expressions and to characterize the 'shape' of the association with different scores (non-monotonicity, non-linearity, complexity ...)}
#'     \item{main assumption}{Two similar scatterplots (btw 2 genes) must have similar "characterization" scores}
#'     \item{main principle}{an association/gene_pair can be represented as a point in a n-dimensional space where n is the number of scores}
#' }
#' @docType package
#' @name mage
NULL
#> NULL