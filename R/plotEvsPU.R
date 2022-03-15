#' Plot error distribution or extremal errors vs prediction uncertainty
#'
#' @param x     (vector) prediction uncertainty
#' @param y     (vector) error
#' @param type  (string) plot type. Default 'points'
#' @param xlab  (string) x axis label. Default 'Pred. unc.'
#' @param xlim  (vector) limits of the x axis
#' @param logX  (logical) log-transform x
#' @param ylab  (string) y axis label. Default 'Error'
#' @param ylim  (vector) limits of the y axis
#' @param nBin  (integer) number of intervals for local coverage stats
#' @param slide (logical) use sliding window for subsetting (x,y)
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param scalePoints (numeric) scale factor for points size
#'
#' @return Plots y vs x or extremal values of abs(y) vs x.
#' @export
#'
#' @examples
plotEvsPU =  function(
  x,
  y,
  type      = c('points','ext'),
  xlab      = 'Error uncertainty',
  ylab      = 'Error',
  logX      = FALSE,
  title     = NULL,
  colPoints = NULL,
  xlim      = NULL,
  ylim      = NULL,
  scalePoints = 0.75,
  nBin      = NULL,
  slide     = TRUE,
  label     = 0,
  gPars
) {

  if (length(x)*length(y) == 0)
    return()

  type = match.arg(type)

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
    yaxs = 'i'
  )

  log = ''
  if(logX)
    log = 'x'
  if(is.null(xlim))
    xlim = range(x)

  if(type == 'points') {
    if(is.null(ylim))
      ylim = 2*max(x)*c(-1,1)

    pch = 16
    colp = cols_tr2[5]

    plot(
      x, y,
      pch = pch,
      col = if(!is.null(colPoints)) colPoints else colp,
      xlim = xlim,
      ylim = ylim,
      xlab = xlab,
      ylab = ylab,
      log = log,
      cex = scalePoints*cex,
      main = title
    )
    grid()
    abline(h = 0, lty = 3)
    for (s in c(-3,-2, -1, 1, 2, 3))
      abline(
        a = 0,
        b = s,
        untf = TRUE,
        lty = 2,
        lwd = 1.5 * lwd,
        col = cols[3]
      )
    box()

  } else {
    if(is.null(nBin))
      nBin  = max(min(floor(N/100),100),10)
    if(nBin <= 0)
      stop('>>> nBin should be > 0')
    if(is.null(slide))
      slide = nBin <= 4

    ord  = order(x)
    xOrd = x[ord]
    yOrd = y[ord]

    # Design local areas
    intrv = ErrViewLib::genIntervals(N, nBin, slide)

    # Local stats
    yext = mint = c()
    for (i in 1:intrv$nbr) {
      sel = intrv$lwindx[i]:intrv$upindx[i]
      # Take max within the interval
      yext[i] = max(abs(yOrd[sel]))
      # Center of the interval
      mint[i] = 0.5*sum(range(xOrd[sel]))
    }

    if(is.null(ylim))
      ylim = c(0,1.05*max(yext))
    if(!any(is.finite(ylim)))
      ylim = xlim

    colp = cols[2]

    plot(
      mint, yext,
      type = 'l', lwd = 2*lwd,
      col = if(!is.null(colPoints)) colPoints else colp,
      xlim = xlim,
      ylim = ylim,
      xlab = xlab,
      ylab = ylab,
      log = log,
      cex = scalePoints*cex,
      main = title
    )
    grid()
    for (s in 1:4)
      abline(
        a = 0,
        b = s,
        untf = TRUE,
        lty = 2,
        lwd = 1.5 * lwd,
        col = cols[3]
      )
    box()
  }

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)

}