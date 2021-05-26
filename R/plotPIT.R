#' Plot PIT histograms
#'
#' @param R (vector) a set of reference values
#' @param C  (vector) a set of predicted values
#' @param uC  (vector) a set of prediction uncertainties
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
plotPIT = function(
  R, C, uC,
  col = 5,
  breaks = 10,
  title = '',
  label = 0,
  gPars
) {

  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = 's',
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  # Uncertainty on uniform histogram
  v  = length(C)/breaks
  uv = sqrt(v)

  h = hist(pnorm((R-C)/uC), breaks = breaks, plot=FALSE)$counts
  hist(
    pnorm((R-C)/uC),
    col=cols_tr2[col],
    xlim = c(0,1),
    xlab = 'PIT',
    ylim = c(0,max(c(v+2*uv,h))),
    breaks = breaks,
    main = title
  )

  abline(h=c(v-2*uv,v,v+2*uv), lty=c(2,1,2), lwd = lwd, col= cols[2])

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)

}