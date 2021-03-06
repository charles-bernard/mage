% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{get_relevant_k}
\alias{get_relevant_k}
\title{get_relevant_k}
\usage{
get_relevant_k(scores_tab, ix = NULL, kmax = 12)
}
\arguments{
\item{scores_tab}{data.table, the table of scores}

\item{ix}{boolean vector of indexes to which \code{scores_tab} must be reduced (sampling output)}

\item{kmax}{maximun number of clusters allowed}
}
\value{
a list of 2 elements
\describe{
  \item{NbClust_output}{A vector of integer storing the best k returned by each method}
  \item{majority_rule_k}{k most often returned}
}
}
\description{
calls NbClust package to return a list of relevant k
according to different methodologies
}
