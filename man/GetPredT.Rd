% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor.R
\name{GetPredT}
\alias{GetPredT}
\title{GetPredT}
\usage{
GetPredT(
  object,
  chains = c(1:10),
  filter_chain = TRUE,
  mean_th = 0.45,
  adjust = TRUE
)
}
\arguments{
\item{object}{mgpfact object}

\item{chains}{The Markov chains for parameter statistics;
if not specified, statistics will be conducted for all}

\item{filter_chain}{Logical value, indicating whether to discard some low-quality chains}

\item{mean_th}{Filtering Markov chains based on the threshold set
for the correlation of pseudotime between different Markov chains}

\item{adjust}{Logical value,
indicating whether to automatically adjust the representation of
differentiation direction by pseudotime}
}
\description{
use this function, obtain correlations under different
iterations of different chains
}
