#' Plot error E vs predicted value V or prediction uncertainty uE
#'
#' @param X (vector) prediction uncertainty, uE, or predicted value, V
#' @param Y (vector) error or z-score
#' @param type (string) type of guide lines ('prop' or 'horiz')
#' @param runExt (logical) plot extremal points
#' @param runQuant (logical) plot running quantiles (cf. probs)
#' @param runMean (logical) plot running mean
#' @param runVar (logical) plot running variance
#' @param runMS (logical) plot running mean squares (exclusive with `runVar`)
#' @param runMode (logical) plot running mode of distribution
#' @param probs (vector) probability levels for quantile-based confidence intervals
#' @param cumMAE (logical) plot cumulative MAE
#' @param xlab (string) x axis label
#' @param xlim (vector) limits of the X axis
#' @param logX (logical) log-transform X
#' @param logY (logical) log-transform Y after taking absolute value
#' @param ylab (string) y axis label
#' @param ylim (vector) limits of the y axis
#' @param nBin (integer) number of intervals for local coverage stats
#' @param slide (logical) use sliding window for subsetting (X,Y)
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param scalePoints (numeric) scale factor for points size
#' @param showLegend (logical) display legend
#' @param legLoc (string) location of legend (see \link[grDevices]{xy.coord})
#'
#' @return Plots E vs uE or V
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)  # Generate errors
#'   plotEvsPU(uE, E, runQuant = TRUE)
#' }
plotEvsPU =  function(
  X, Y,
  type      = c('prop','horiz'),
  runExt    = FALSE,
  runQuant  = FALSE,
  runMode   = FALSE,
  runMean   = FALSE,
  runVar    = FALSE,
  runMS     = FALSE,
  cumMAE    = FALSE,
  probs     = c(0.95),
  xlab      = NULL,
  ylab      = NULL,
  logX      = FALSE,
  logY      = FALSE,
  title     = NULL,
  myColors  = c(5,4,3,1,2),
  xlim      = NULL,
  ylim      = NULL,
  scalePoints = 0.5,
  nBin      = NULL,
  slide     = TRUE,
  label     = 0,
  showLegend = TRUE,
  legLoc    = 'topleft',
  gPars     = ErrViewLib::setgPars()
) {

  getmode = function(X) {
    d <- density(X)
    d$x[which.max(d$y)]
  }

  if (as.numeric(length(X))*as.numeric(length(Y)) == 0)
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

  if(logY)
    Y = log10(abs(Y))

  if(runExt | runQuant | runMode | runVar | runMean) {
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
    ymin = ymax = ymode = yvar = ymean = mint = c()
    nq = length(probs)
    qmin = qmax = matrix(NA,ncol=nq,nrow=intrv$nbr)
    for (i in 1:intrv$nbr) {

      # Subsetting
      sel = intrv$lwindx[i]:intrv$upindx[i]
      x = xOrd[sel]
      y = yOrd[sel]

      # Take min & max within the interval
      ymin[i] = min(y)
      ymax[i] = max(y)
      if(runQuant) {
        for (j in seq_along(probs)) {
          p = probs[j]
          alpha = c((1 - p)/2, (1 + p)/2)
          q = ErrViewLib::vhd(y, p = alpha)
          qmin[i, j] = q[1]
          qmax[i, j] = q[2]
        }
      }

      if(runMode)
        ymode[i] = getmode(y)

      if(runVar)
        yvar[i] = var(y)

      if(runMS)
        yvar[i] = mean(y^2)

      if(runMean)
        ymean[i] = mean(y)

      # Center of the interval
      mint[i] = 0.5*sum(range(x))
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

  if(is.null(xlim)){
    if(type == 'prop')
      xlim = c(max(0,min(X)), 1.1*max(X))
    else
      xlim = range(X)
  }

  if(is.null(ylim)){
    if(logY)
      ylim = range(Y)
    else
      ylim = c(-max(abs(Y)), max(abs(Y)))
  }

  if(is.null(xlab))
    xlab = ifelse(
      type == 'horiz',
      'Predicted Value, V',
      'Prediction Uncertainty, uE'
    )

  if(is.null(ylab))
    ylab = ifelse(
      type == 'horiz',
      'Z-score',
      'Error, E'
    )

  pch = 16

  colp = cols_tr2[myColors[1]] # Points
  cole = cols[myColors[2]]     # Extremal lines
  colr = cols[myColors[3]]     # Quantile reg
  colg = cols[myColors[4]]     # Guide lines
  colc = cols[myColors[5]]     # Cum. MAE
  colm = cols[myColors[3]]     # Mode
  colv = cols[myColors[5]]     # Variance
  cola = cols[myColors[3]]     # Variance

  plot(
    X, Y,
    type = 'n',
    xlim = xlim,
    xaxs = 'i',
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    log = log,
    main = title
  )
  grid()
  # Guide lines
  if(logY) {
    lines(
      xOrd,
      log10(xOrd),
      lty = 2,
      col = colg
    )
  } else {
    abline(h = 0, lty = 3)
    klev = c(-3,-2, -1, 1, 2, 3)
    for (s in klev)
      if(type == 'horiz') {
        abline(
          h = s,
          lty = 2,
          col = colg
        )
      } else {
        abline(
          a = 0,
          b = s,
          untf = TRUE,
          lty = 2,
          col = colg
        )
      }
  }
  points(
    X, Y,
    pch = pch,
    col = colp,
    cex = scalePoints*cex
  )

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
    for (j in seq_along(probs)){
      matlines(
        mint, cbind(qmin[,j],qmax[,j]),
        lty = 1,
        lwd = 2*lwd,
        col = colr)
    }
    klev = c(-2, 2)
    for (s in klev)
      if(type == 'horiz') {
        abline(
          h = s,
          lty = 2,
          lwd = 2*lwd,
          col = colr
        )
      } else {
        abline(
          a = 0,
          b = s,
          untf = TRUE,
          lty = 2,
          lwd = 2*lwd,
          col = colr
        )
      }
    legs = c(legs,'Quantiles')
    pleg = c(pleg,NA)
    cleg = c(cleg,colr)
    tleg = c(tleg,1)
  }

  if(runMode) {
    lines(
      mint, ymode,
      lty = 1,
      lwd = 2*lwd,
      col = colm)
    legs = c(legs,'Mode')
    pleg = c(pleg,NA)
    cleg = c(cleg,colm)
    tleg = c(tleg,1)
  }

  if(runMean) {
    lines(
      mint, ymean,
      lty = 1,
      lwd = 2*lwd,
      col = cola)
    legs = c(legs,'Mean')
    pleg = c(pleg,NA)
    cleg = c(cleg,cola)
    tleg = c(tleg,1)
  }

  if(runVar | runMS) {
    lines(
      mint, yvar,
      lty = 1,
      lwd = 2*lwd,
      col = colv)
    klev = c(1)
    for (s in klev)
      if(type == 'horiz') {
        abline(
          h = s,
          lty = 2,
          lwd = 2*lwd,
          col = colv
        )
      } else {
        abline(
          a = 0,
          b = s,
          untf = TRUE,
          lty = 2,
          lwd = 2*lwd,
          col = colv
        )
      }
    legs = c(legs,ifelse(runVar,'Variance','Mean squares'))
    pleg = c(pleg,NA)
    cleg = c(cleg,colv)
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

  if (showLegend)
    legend(
      legLoc, bg = 'white', box.col = 'white',
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