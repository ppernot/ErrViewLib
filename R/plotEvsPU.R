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
#' @param slide (logical) use sliding window for subsetting (X,Z)
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
  xlab      = 'Pred. unc.',
  ylab      = 'Error',
  logX      = FALSE,
  title     = NULL,
  colPoints = NULL,
  xlim      = NULL,
  ylim      = NULL,
  scalePoints = 0.75,
  nBin      = NULL,
  slide     = NULL,
  gPars
) {

  if (length(x)*length(y) == 0)
    return()

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mgp = mgp,
    pty = 's',
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  pch = 16
  colp = cols_tr2[5]

  log = ''
  if(logX)
    log = 'x'
  if(is.null(xlim))
    xlim = range(x)

  if(type == 'points') {
    if(is.null(ylim))
      ylim = 2*max(x)*c(-1,1)

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
    suC = sort(x)
    lines(
      suC, 2*suC,
      lty = 2,
      lwd = lwd,
      col = cols[3]
    )
    lines(
      suC, -2*suC,
      lty = 2,
      lwd = lwd,
      col = cols[3]
    )
    lines(
      suC, suC,
      lty = 2,
      lwd = lwd,
      col = cols[3]
    )
    lines(
      suC, -suC,
      lty = 2,
      lwd = lwd,
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
    if (slide) {
      # Sliding interval of width nLoc
      nLoc = floor(N / nBin)
      # Nbr of intervals
      nbr  = N - nLoc +1
      # Lower index of interval in ordered data
      lwindx = 1:nbr
      # Upper index
      upindx = lwindx + nLoc -1

    } else {
      # Breakpoints of nearly equi-sized C intervals
      p    = seq(0, 1, length.out = nBin + 1)[1:nBin]
      br   = ErrViewLib::vhd(xOrd, p = p)
      # Nbr of intervals
      nbr  = length(br)
      # Lower index of interval in ordered data
      lwindx = upindx = c()
      lwindx[1] = 1
      for (i in 2:nbr)
        lwindx[i] = which(xOrd > br[i])[1]
      # Upper index
      for (i in 1:(nbr-1))
        upindx[i] = lwindx[i+1]-1
      upindx[nbr] = N

      if(min(upindx-lwindx) < N/nBin/2 |
         sum(upindx-lwindx+1) != N      )
        stop('>>> Pb in equi-sized intervals design')
    }
    yext = mint = c()
    for (i in 1:nbr) {
      sel = lwindx[i]:upindx[i]
      M = length(sel)
      yLoc = yOrd[sel]
      yext[i] = max(abs(yLoc))
      mint[i] = 0.5*sum(range(xOrd[sel]))
    }

    if(is.null(ylim))
      ylim = c(0,max(yext))

    plot(
      mint, yext,
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
    suC = sort(x)
    lines(
      suC, 3*suC,
      lty = 2,
      lwd = lwd,
      col = cols[3]
    )
    lines(
      suC, 2*suC,
      lty = 2,
      lwd = lwd,
      col = cols[3]
    )
    lines(
      suC, suC,
      lty = 2,
      lwd = lwd,
      col = cols[3]
    )
    box()
  }
}