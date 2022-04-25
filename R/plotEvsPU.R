#' Plot error E vs predicted value V or prediction uncertainty uE
#'
#' @param X (vector) prediction uncertainty, uE, or predicted value, V
#' @param Y (vector) error or z-score
#' @param type (string) type of guide lines ('prop' or 'horiz')
#' @param runExt (logical) plot extremal points
#' @param runQuant (logical) plot quantile regression
#' @param cumMAE (logical) plot cumulative MAE
#' @param xlab (string) x axis label
#' @param xlim (vector) limits of the x axis
#' @param logX (logical) log-transform x
#' @param ylab (string) y axis label
#' @param ylim (vector) limits of the y axis
#' @param nBin (integer) number of intervals for local coverage stats
#' @param slide (logical) use sliding window for subsetting (X,Y)
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param scalePoints (numeric) scale factor for points size
#'
#' @return Plots E vs uE or V
#' @export
#'
#' @examples
plotEvsPU =  function(
  X, Y,
  type      = c('prop','horiz'),
  runExt    = FALSE,
  runQuant  = FALSE,
  cumMAE    = FALSE,
  xlab      = NULL,
  ylab      = NULL,
  logX      = FALSE,
  title     = NULL,
  myColors  = c(5,4,3,2,7),
  xlim      = NULL,
  ylim      = NULL,
  scalePoints = 0.5,
  nBin      = NULL,
  slide     = TRUE,
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  if (length(X)*length(Y) == 0)
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
    yaxs = 'i',
    cex.main = 1
  )

  if(runExt | runQuant) {
    N = length(Y)

    if(is.null(nBin))
      nBin  = ceiling(2*N^0.333)

    if(nBin <= 0)
      stop('>>> nBin should be > 0')

    ord  = order(X)
    xOrd = X[ord]
    yOrd = Y[ord]

    # Design local areas
    intrv = ErrViewLib::genIntervals(N, nBin, slide)

    # Local stats
    ymin = ymax = qmin = qmax = mint = c()
    for (i in 1:intrv$nbr) {
      # Subsetting
      sel = intrv$lwindx[i]:intrv$upindx[i]
      x = xOrd[sel]
      y = yOrd[sel]
      # Take min & max within the interval
      ymin[i] = min(y)
      ymax[i] = max(y)
      if(runQuant) {
        q = ErrViewLib::vhd(y)
        qmin[i] = q[1]
        qmax[i] = q[2]
      }
      # Center of the interval
      mint[i] = 0.5*sum(range(xOrd[sel]))
    }
  }

  if(cumMAE) {
    ord  = order(X)
    xOrd = X[ord]
    yOrd = Y[ord]
    cMAE = cumsum(abs(yOrd)) / 1:length(yOrd)
  }

  log = ''
  if(logX)
    log = 'x'

  if(is.null(xlim))
    xlim = range(X)

  if(is.null(ylim))
    ylim = c(-max(abs(Y)), max(abs(Y)))

  if(is.null(xlab))
    xlab = ifelse(
      type == 'PV',
      'Predicted Value, V',
      'Prediction Uncertainty, uE'
    )

  if(is.null(ylab))
    ylab = ifelse(
      type == 'PV',
      'z-score, E/uE',
      'Error, E'
    )

  pch = 16

  colp = cols_tr2[myColors[1]] # Points
  cole = cols[myColors[2]]     # Extremal lines
  colr = cols[myColors[3]]     # Quantile reg
  colg = cols[myColors[4]]     # Guide lines
  colc = cols[myColors[5]]     # Cum. MAE

  plot(
    X, Y,
    pch = pch,
    col = colp,
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    log = log,
    cex = scalePoints*cex,
    main = title
  )
  # grid()
  # Guide lines
  abline(h = 0, lty = 3)
  for (s in c(-3,-2, -1, 1, 2, 3))
    if(type == 'horiz') {
      abline(
        h = s,
        lty = 2,
        lwd = 2*lwd,
        col = colg
      )
    } else {
      abline(
        a = 0,
        b = s,
        untf = TRUE,
        lty = 2,
        lwd = 2*lwd,
        col = colg
      )
    }

  legs = 'Data'
  pleg = pch
  cleg = colp
  tleg = NA

  if(runExt) {
    matlines(
      mint, cbind(ymin,ymax),
      lty = 1,
      lwd = 2*lwd,
      col = cole)
    legs = c(legs,'Extrema')
    pleg = c(pleg,NA)
    cleg = c(cleg,cole)
    tleg = c(tleg,1)
  }

  if(runQuant) {
    matlines(
      mint, cbind(qmin,qmax),
      lty = 1,
      lwd = 2*lwd,
      col = colr)
    legs = c(legs,'Quantiles')
    pleg = c(pleg,NA)
    cleg = c(cleg,colr)
    tleg = c(tleg,1)
  }

  if(cumMAE) {
    matlines(
      xOrd,cMAE,
      lty = 1,
      lwd = 2*lwd,
      col = colc)
    legs = c(legs,'Cum. MAE')
    pleg = c(pleg,NA)
    cleg = c(cleg,colc)
    tleg = c(tleg,1)
  }

  legend(
    'topleft', bg = 'white', box.col = 'white',
    legend = legs,
    pch = pleg,
    lty = tleg,
    col = cleg,
    lwd = 2*lwd,
    cex = 0.9
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