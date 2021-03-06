% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operations_on_scores.R
\name{standardize_scores}
\alias{standardize_scores}
\title{standardize_scores}
\usage{
standardize_scores(x)
}
\arguments{
\item{x}{data.frame, table of scores}
}
\value{
returns a standardized table of scores
}
\description{
Standardize each column of the table of scores via the Zscore: \deqn{x = ( x - mean(column) ) / ( std(column) )}
}
