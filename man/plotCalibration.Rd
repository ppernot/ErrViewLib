% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCalCurve.R
\name{plotCalibration}
\alias{plotCalibration}
\title{Plot calibration curve}
\usage{
plotCalibration(
  E,
  uE,
  prob = seq(0.005, 0.995, by = 0.005),
  dist = c("Normal", "Student", "Uniform", "Laplace", "Normp", "T4", "Normp4"),
  shape = 4,
  score = TRUE,
  plot = TRUE,
  title = "",
  label = 0,
  gPars = ErrViewLib::setgPars()
)
}
\arguments{
\item{E}{(vector) a vector of errors}

\item{uE}{(vector) a vector of uncertainties}

\item{dist}{(string) a distribution (default: `Normal`)}

\item{shape}{(numeric) shape parameter for the T and Normp distributions}

\item{score}{(logical) evaluate the calibration scores (default: `TRUE`)}

\item{plot}{(logical) plot the results (default: `TRUE`)}

\item{title}{(string) a title to display above the plot}

\item{label}{(integer) index of letter for subplot tag}

\item{gPars}{(list) graphical parameters}
}
\value{
Invisibly returns a list of calibration statistics. Mainly used
  for its plotting side effect.
}
\description{
Plot calibration curve
}
