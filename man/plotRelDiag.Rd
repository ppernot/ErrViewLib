% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotRelDiag.R
\name{plotRelDiag}
\alias{plotRelDiag}
\title{Plot reliability diagram}
\usage{
plotRelDiag(
  X,
  Y,
  ordX = NULL,
  aux = NULL,
  logX = FALSE,
  nBin = NULL,
  equiPop = TRUE,
  popMin = 30,
  logBin = TRUE,
  intrv = NULL,
  method = c("bootstrap", "cho"),
  BSmethod = c("bca", "perc", "basic"),
  nBoot = 500,
  slide = FALSE,
  plot = TRUE,
  score = TRUE,
  xlim = NULL,
  ylim = NULL,
  title = "",
  unit = "",
  label = 0,
  add = FALSE,
  col = 5,
  gPars = ErrViewLib::setgPars()
)
}
\arguments{
\item{X}{(vector) vector of uncertainties}

\item{Y}{(vector) vector of errors}

\item{ordX}{(vector) set of abscissas to order sample}

\item{aux}{(vector) auxilliary vector to resolve ties in sample sorting}

\item{logX}{(logical) log-transform X}

\item{nBin}{(integer) number of intervals for local coverage stats}

\item{equiPop}{(logical) generate intervals  with equal bin counts
(default: `equiPop = TRUE`)}

\item{popMin}{(integer) minimal bin count in an interval}

\item{logBin}{(logical) if `equiPop = FALSE`, one can choose between
equal range intervals, or equal log-range intervals
(default `logBin = TRUE`)}

\item{intrv}{(object) intervals generated by `genIntervals` (default: `NULL`)}

\item{method}{(string) one of 'bootstrap' (default) and 'cho'}

\item{BSmethod}{(string) bootstrap variant}

\item{nBoot}{(integer) number of bootstrap replicas}

\item{slide}{(logical) use sliding window for subsetting (X,Z)}

\item{plot}{(function) produce plot ?}

\item{score}{(logical) estimate calibration metrics ?}

\item{xlim}{(vector) min and max values of X axis}

\item{ylim}{(vector) limits of the y axis}

\item{title}{(string) a title to display above the plot}

\item{unit}{(string) unit string to add to xlab and ylab}

\item{label}{(integer) index of letter for subplot tag}

\item{add}{(logical) add to previous graph ?}

\item{col}{(integer) color index of curve to add}

\item{gPars}{(list) graphical parameters}
}
\value{
Used for its plotting side effect.
}
\description{
Plot reliability diagram
}
\examples{
\donttest{
  uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
  E   = rnorm(uE, mean=0, sd=uE)  # Generate errors
  plotRelDiag(uE, E, nBin = 6, nBoot = 500)
}
}
