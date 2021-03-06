% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDeltaCDF.R
\name{plotDeltaCDF}
\alias{plotDeltaCDF}
\title{Plot of the CDF of the differences of a pair of absolute errors samples}
\usage{
plotDeltaCDF(
  err,
  meth1,
  meth2,
  eps = NULL,
  xmin = NULL,
  xmax = NULL,
  xlab = NULL,
  units = "a.u.",
  main = "",
  nboot = 1000,
  label = 0,
  showSIP = TRUE,
  showMLG = TRUE,
  showDmue = TRUE,
  showCI = TRUE,
  gPars
)
}
\arguments{
\item{gPars}{}
}
\value{

}
\description{
Plots also auxiliary statistics (SIP, MG, ML) with bootstrapped CIs.
}
