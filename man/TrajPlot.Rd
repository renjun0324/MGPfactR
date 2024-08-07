% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/track.R
\name{TrajPlot}
\alias{TrajPlot}
\title{TrajPlot}
\usage{
TrajPlot(
  object,
  box_fill = NULL,
  plot_title = TRUE,
  g_title = NULL,
  col = NULL,
  plot_tb = TRUE,
  tb_pred = NULL,
  pointSize = 5,
  pointAlpha = 0.4,
  pointLabel = FALSE,
  pointLabelsize = 3,
  rug = TRUE,
  rugAlpha = 0.3,
  legend = TRUE,
  legend_title = NULL,
  lm = FALSE,
  lineMethod = loess,
  lineType = "solid",
  lineSize = 1,
  lineAlpha = 0.3,
  se = F,
  span = 0.75,
  formula = y ~ x,
  family = NULL,
  vlineType = "dashed",
  ...
)
}
\arguments{
\item{object}{MGPfact object}

\item{box_fill}{box color}

\item{col}{the colors of different branches}

\item{tb_pred}{bifurcation point}

\item{pointSize}{point size}

\item{pointAlpha}{point alpha}

\item{pointLabel}{logical value, whether label point}

\item{pointLabelsize}{label size}

\item{rug}{logical value, whether to plot rug}

\item{rugAlpha}{rug alpha}

\item{legend}{locical value, whether to add legend}

\item{legend_title}{legend title}

\item{lm}{whether add lm line}

\item{lineMethod}{linear model methods}

\item{lineType}{line type}

\item{lineSize}{line size}

\item{lineAlpha}{line alpha}

\item{save}{logical value, whether to save pdf}

\item{title}{plot title}

\item{tb}{logical value, whether label bifurcation line}
}
\description{
smooth trajectory visualization
}
