% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotRankMat.R
\name{plotRankMat}
\alias{plotRankMat}
\title{Plots of Ranking Probability matrix Pr}
\usage{
plotRankMat(
  E,
  tab = NULL,
  score = "mue",
  type = "levels",
  method = "square",
  nMC = 1000,
  cex.lab = 1,
  show.main = TRUE,
  offset = 0.7,
  M = nrow(E),
  label = 0,
  gPars = ErrViewLib::setgPars()
)
}
\arguments{
\item{E}{-}

\item{tab}{-}

\item{score}{-}

\item{type}{-}

\item{method}{-}

\item{nMC}{-}

\item{cex.lab}{-}

\item{show.main}{-}

\item{offset}{-}

\item{M}{-}

\item{gPars}{-}
}
\description{
Interface to \code{corrplot::corrplot()} or mode/CI graph
}
