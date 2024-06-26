% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotLZV.R
\name{plotLZMS}
\alias{plotLZMS}
\title{Plot local Mean Squared z-score <Z^2> to assess calibration and tightness}
\usage{
plotLZMS(
  X,
  Z,
  aux = NULL,
  varZ = 1,
  logX = FALSE,
  nBin = NULL,
  equiPop = TRUE,
  popMin = 100,
  logBin = TRUE,
  intrv = NULL,
  plot = TRUE,
  slide = FALSE,
  nBoot = 5000,
  parallel = FALSE,
  method = c("bootstrap", "stud", "auto"),
  BSmethod = c("bca", "perc", "basic", "stud"),
  xlab = "Conditioning variable",
  skewup = 0.8,
  xlim = NULL,
  ylim = NULL,
  title = "",
  score = FALSE,
  add = FALSE,
  col = 5,
  colInv = 2,
  label = 0,
  gPars = ErrViewLib::setgPars()
)
}
\arguments{
\item{X}{(vector) abscissae of the Z values}

\item{Z}{(vector) set of z-score values to be tested}

\item{aux}{(vector) auxilliary vector to resolve ties in X sorting}

\item{varZ}{(numeric) target value for Var(Z) (default `1`)}

\item{logX}{(logical) log-transform X}

\item{nBin}{(integer) number of intervals for local coverage stats}

\item{equiPop}{(logical) generate intervals  with equal bin counts
(default: `equiPop = TRUE`)}

\item{popMin}{(integer) minimal bin count in an interval}

\item{logBin}{(logical) if `equiPop = FALSE`, one can choose between
equal range intervals, or equal log-range intervals
(default `logBin = TRUE`)}

\item{intrv}{(object) intervals generated by `genIntervals` (default: `NULL`)}

\item{plot}{(logical) plot the results}

\item{slide}{(logical) use sliding window for subsetting (X,Z)}

\item{nBoot}{(integer) number of bootstrap replicas}

\item{parallel}{(logical) parallelized bootstrap (default: `FALSE`)}

\item{method}{(string) method used to estimate 95 percent CI on <Z^2>}

\item{BSmethod}{(string) bootstrap variant}

\item{xlab}{(string) X axis label}

\item{skewup}{(numeric) upper limit for robust skewness of Z^2, used
for reliability estimation (defaul: 0.8). The unreliable results
are grayed out. Set to NULL to inactivate this check.}

\item{xlim}{(vector) min and max values of X axis}

\item{ylim}{(vector) limits of the y axis}

\item{title}{(string) a title to display above the plot}

\item{score}{(logical) estimate calibration stats (default: `FALSE`)}

\item{add}{(logical) add to previous graph ?}

\item{col}{(integer) color index of bin stats}

\item{colInv}{(integer) color index of invalid bin stats}

\item{label}{(integer) index of letter for subplot tag}

\item{gPars}{(list) graphical parameters}
}
\value{
Invisibly returns a list of LZMS results. Mainly used
  for its plotting side effect.
}
\description{
Plot local Mean Squared z-score <Z^2> to assess calibration and tightness
}
\examples{
\donttest{
  uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
  E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
  plotLZMS(uE, E/uE, method = 'cho', ylim = c(0,2))
}
}
