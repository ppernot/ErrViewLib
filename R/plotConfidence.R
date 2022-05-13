#' Plot confidence curve for (uE,E) set
#'
#' @param E (vector) prediction uncertainty, uE, or predicted value, V
#' @param uE (vector) error or z-score
#' @param stat (function) statistic to use
#' @param oracle (logical) plot oracle curve
#' @param col (integer) color index for the curve
#' @param add (logical) add confidence curve to previous plot
#' @param xlab (string) x axis label
#' @param xlim (vector) limits of the x axis
#' @param ylab (string) y axis label
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#'
#' @return Plot confidence curve for (uE,E) set

#' @export
#'
plotConfidence = function(
  E, uE,
  stat   = ErrViewLib::mue,
  oracle = TRUE,
  col    = 2,
  add    = FALSE,
  xlab   = '% discarded',
  xlim   = NULL,
  ylab   = 'MUE',
  ylim   = NULL,
  title  = NULL,
  label  = 0,
  gPars  = ErrViewLib::setgPars()
) {

  if (as.numeric(length(E))*as.numeric(length(uE)) == 0)
    return()

  # Expose gPars list
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
    yaxs = 'i',
    cex.main = 1
  )

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

  vstat = vora = rep(NA,length(pcVec))
  for (i in seq_along(pcVec)) {
    k = pcVec[i]
    sel = 1:floor(k * M / 100)
    if(length(sel) == 0)
      break()
    vstat[i] = stat(E[-sel])
    if(oracle)
      vora[i]  = stat(O[-sel])
  }
  # Add 0% point
  pcVec = c(0,pcVec)
  vstat = c(stat(E),vstat)
  if(oracle)
    vora = c(stat(O),vora)

  if(add) {
    lines(pcVec, vstat, col = cols[col])

  } else {

    if (is.null(xlim))
      xlim = range(pcVec)

    if (is.null(ylim))
      ylim = c(0, 1.1 * max(vstat, na.rm = TRUE))

    plot(
      pcVec, vstat,
      type = 'l',
      col  = cols[col],
      xlab = xlab,
      xlim = xlim,
      ylab = ylab,
      ylim = ylim,
      main = title
    )
    grid()
    if(oracle)
      lines(pcVec, vora, lty=2, col=cols[1])

    legend(
      'bottomleft', bty = 'n', inset = 0.05,
      legend = c('Oracle','Data'),
      lty = c(2,1),
      col = cols[c(1,col)],
      pch = NA
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