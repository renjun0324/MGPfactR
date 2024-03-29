% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{MURPDownsampling}
\alias{MURPDownsampling}
\title{MURPDownsampling}
\usage{
MURPDownsampling(
  object,
  omega = 0.5,
  iter = 10,
  fast = T,
  max_murp = 500,
  pca.center = FALSE,
  pca.scale = FALSE,
  cores = 1,
  seed = 723,
  plot = T
)
}
\arguments{
\item{object}{MGPfact object}

\item{omega}{MURP parameter: omega value in pseudo-BIC calculating}

\item{iter}{MURP parameter: iter the number of iterations}

\item{fast}{MURP parameter:speed up calculations}

\item{max_murp}{MURP parameter: max downsampling number}

\item{pca.center}{MURPPCA parameter: a logical value indicating whether the variables should be
shifted to be zero centered}

\item{pca.scale}{MURPPCA parameter: a logical value indicating whether the variables should be
scaled to have unit variance before the analysis takes place}

\item{cores}{the number of cores to use in multi-threads calculating, default is 1.
Multi-threaded calculation is not recommended except for users with server support}

\item{seed}{random seed, default is 723}

\item{plot}{PCANestedGridPlot}
}
\description{
MURP downsampling function
}
