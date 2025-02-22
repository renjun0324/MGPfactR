% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimensional_reduction.R
\name{RUNPCA}
\alias{RUNPCA}
\title{RUNPCA}
\usage{
RUNPCA(object, pc_num = 30, center = TRUE, scale = FALSE)
}
\arguments{
\item{object}{MGPfact object}

\item{pc_num}{Number of principal components}

\item{center}{a logical value indicating whether the variables should be shifted to be zero centered.
Alternately, a vector of length equal the number of columns of x can be supplied. The value is passed to scale.}

\item{scale}{a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.}
}
\description{
Perform principal component analysis on the original gene expression matrix
}
