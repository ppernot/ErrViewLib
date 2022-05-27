#' Plot PIT histograms
#'
#' @param Z (vector) a set of z-zscores
#' @param dist (string) a distribution
#' @param shape (numeric) shape parameter (> 0, maybe non-integer).
#' @param col (integer) a color number for the histigram
#' @param breaks (integer) number of breaks in the histogram
#' @param title (string) a title for the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#'
#' @return a plot
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   plotPIT(E/uE)
#' }
plotPIT = function(
  Z,
  dist  = c('norm','t'),
  shape = 2,
  col = 5,
  breaks = 10,
  title = '',
  label = 0,
  gPars = ErrViewLib::setgPars()
) {

  dist = match.arg(dist)

  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = 's',
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    xaxs = 'i',
    yaxs = 'i'
  )

  # CI on uniform histogram
  v  = length(Z)/breaks
  lwr = qpois(0.025, v)
  upr = qpois(0.975, v)

  # Probability Integral Transform using unit variance hyp.
  pit = switch(
    dist,
    norm = normalp::pnormp(
      Z * sqrt(shape^(2/shape)*gamma(3/shape)/gamma(1/shape)),
      p = shape) ,
    t    = pt(
      Z * sqrt(shape / (shape-2)),
      df = shape)
  )

  # Histogram
  h = hist(pit, breaks = breaks, plot=FALSE)$counts # To estimate ylim
  hist(
    pit,
    col=cols_tr[col],
    xlim = c(0,1),
    xlab = 'PIT',
    ylim = c(0,1.1*max(c(upr,h))),
    breaks = breaks,
    main = title
  )

  # 95% CI
  abline(
    h=c(lwr,v,upr),
    lty=c(2,1,2),
    lwd = lwd,
    col= cols[2]
  )

  box()

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)

}