#' Plot confidence curve for (uE,E) set
#'
#' @param E (vector) prediction uncertainty, uE, or predicted value, V
#' @param uE (vector) error or z-score
#' @param stat (function) statistic to use
#' @param oracle (logical) plot oracle curve
#' @param col (integer) color index for the curve
#' @param type (string) curve type: line ('l') or points ('p')
#' @param add (logical) add confidence curve to previous plot
#' @param xlab (string) x axis label
#' @param xlim (vector) limits of the x axis
#' @param ylab (string) y axis label
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param legend (string) legend for the dataset
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#'
#' @return Plot confidence curve for (E,uE) set

#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)  # Generate errors
#'   plotConfidence(E,uE)
#' }
#'
plotConfidence = function(
  E, uE,
  stat   = ErrViewLib::mue,
  oracle = TRUE,
  col    = 2,
  type   = c('l','p'),
  add    = FALSE,
  xlab   = 'k% discarded',
  xlim   = NULL,
  ylab   = 'MAE / MAE0',
  ylim   = NULL,
  title  = NULL,
  label  = 0,
  legend = NULL,
  gPars  = ErrViewLib::setgPars()
) {

  type = match.arg(type)

  if (as.numeric(length(E))*as.numeric(length(uE)) == 0)
    return()

  # Reorder data
  io  = order(uE, decreasing = TRUE)
  uE  = uE[io]
  E   = E[io]
  if(oracle) {
    io = order(abs(E), decreasing = TRUE) # Perfect set for oracle
    O  = E[io]
  }

  M = length(uE)
  pcVec = 1:100 # Vector of percentages

  S0 = stat(E)
  vstat = vora = rep(NA,length(pcVec))
  for (i in seq_along(pcVec)) {
    k = pcVec[i]
    sel = 1:floor(k * M / 100)
    if(length(sel) == 0)
      break()
    vstat[i] = stat(E[-sel]) / S0
    if(oracle)
      vora[i]  = stat(O[-sel]) / S0
  }
  # Add 0% point
  pcVec = c(0,pcVec)
  vstat = c(1,vstat)
  if(oracle)
    vora = c(1,vora)

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  if(add) {

    lines(pcVec, vstat,
          type = type,
          pch = 16,
          lwd = 2 * lwd,
          col = cols[col])

  } else {

    par(
      mfrow = c(1, 1),
      mar = mar,
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd,
      xaxs = 'i',
      yaxs = 'i',
      cex.main = 1
    )

    if (is.null(xlim))
      xlim = range(pcVec)

    if (is.null(ylim))
      ylim = c(0, 1)

    plot(
      pcVec, vstat,
      type = type,
      lty  = 1,
      pch  = 16,
      lwd  = 2*lwd,
      col  = cols[col],
      xlab = xlab,
      xlim = xlim,
      ylab = ylab,
      ylim = ylim,
      main = title
    )
    grid()
    if(oracle)
      lines(pcVec, vora, lty=2, lwd = 2*lwd, col=cols[1])

    lty = if(type == 'l') 1 else NA
    pch = if(type == 'l') NA else 16
    if(oracle) {
      lty = c(2,lty)
      pch = c(NA,pch)
    }
    if(!is.null(legend))
      legend(
        'bottomleft', bty = 'n', inset = 0.05,
        legend = if(oracle) c('Oracle',legend) else legend,
        lty = lty,
        lwd = 2*lwd,
        col = if(oracle) cols[c(1,col)] else cols[col],
        pch = pch
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
}