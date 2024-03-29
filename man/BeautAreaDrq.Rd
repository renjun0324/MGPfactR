% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene.R
\name{BeautAreaDrq}
\alias{BeautAreaDrq}
\title{BeautAreaMulti}
\usage{
BeautAreaDrq(
  df,
  mdf,
  g,
  cut_tb = FALSE,
  cols = c("#374e55ff", "#df8f44ff"),
  y_title_vjust = -55,
  max_str = 50
)
}
\arguments{
\item{df}{Data frame containing all genes and trajectory information,
transformed into a long matrix with separate columns for pseudotime and branch label.}

\item{mdf}{dataframe for any one gene on any one trajectory,
including T, branch label, celltype, gene expression}

\item{g}{any gene name}

\item{cut_tb}{logical value, Whether to remove the section before the bifurcation}

\item{cols}{color of gene expression density in different branches}

\item{y_title_vjust}{Distance of y-axis labels from the panel.}

\item{max_str}{Maximum width for cell type strings}
}
\description{
Display the expression density variation of any one gene across the time domain,
combined with the density variation of cell types
}
