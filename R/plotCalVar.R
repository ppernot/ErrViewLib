#' Plot variance calibration curve
#'
#' @param nBin  (integer) number of intervals for local coverage stats
#' @param ylim  (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param X     (vector) vector of uncertainties
#' @param Y     (vector) vector of errors
#' @param logX  (logical) log-transform X
#' @param slide (logical) use sliding window for subsetting (X,Z)
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
#' @param xlim (vector) min and max values of X axis
#' @param add (logical) add to previous graph ?
#' @param col (interger) color index of curve to add
#'
#' @return Used for its plotting side effect.
#' @export
#'
#' @examples
plotCalVar = function(
  X, Y,
  logX      = FALSE,
  nBin      = NULL,
  slide     = NULL,
  BSmethod  = c('bca','perc','basic'),
  nBoot     = 0,
  xlim      = NULL,
  ylim      = NULL,
  title     = '',
  label     = 0,
  add       = FALSE,
  col       = 5,
  gPars     = ErrViewLib::setgPars()
) {

  BSmethod = match.arg(BSmethod)

  N = length(Y)

  if(is.null(nBin))
    nBin  = max(min(floor(N/150),15),2)
  if(nBin <= 0)
    stop('>>> nBin should be > 0')
  if(is.null(slide))
    slide = nBin <= 4

  ord  = order(X)
  xOrd = X[ord]
  yOrd = Y[ord]

  # Design local areas
  intrv = ErrViewLib::genIntervals(N, nBin, slide)
  xl = log(xOrd)
  step = diff(range(xl))/nBin
  lims = exp(min(xl) + (1:nBin -1)* step)
  lwindx = upindx = c()
  lwindx[1] = 1
  for(i in 2:nBin) {
    lwindx[i] = which(xOrd >= lims[i])[1]
    upindx[i-1] = lwindx[i] - 1
  }
  upindx[nBin] = N
  intrv = list(
    lwindx = lwindx,
    upindx = upindx,
    nbr = nBin
  )

  # Local variance values
  mV = loV = upV = mint = c()
  for (i in 1:intrv$nbr) {
    sel = intrv$lwindx[i]:intrv$upindx[i]
    M = length(sel)
    if (M < 5) {
      mV[i] = NA
      next
    }
    yLoc = yOrd[sel]
    mV[i]  = sd(yLoc)
    if(nBoot != 0) {
      fsd = function(x,w)
        sqrt(Hmisc::wtd.var(x,w,normwt=TRUE))
      bst = boot::boot(yLoc, fsd, R = nBoot, stype = 'f')
      bci = boot::boot.ci(bst, conf = 0.95, type = BSmethod)
      ci = switch(
        BSmethod,
        bca   = bci$bca[1, 4:5],
        perc  = bci$perc[1, 4:5],
        basic = bci$basic[1, 4:5]
      )
      loV[i] = ci[1]
      upV[i] = ci[2]
    } else {
      loV[i] = NA
      upV[i] = NA
    }
    mint[i] = sqrt(mean(xOrd[sel]^2))
  }
  # Remove empty intervals
  sel = !is.na(mV)
  mV = mV[sel]
  loV = loV[sel]
  upV = upV[sel]
  mint = mint[sel]


  if(length(gPars) == 0)
    gPars = ErrViewLib::setgPars()

  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  if(!add) {
    par(
      mfrow = c(1, 1),
      mar = mar,
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd,
      cex.main = 1
    )

    if(is.null(xlim))
      if(!any(is.na(loV)))
        xlim = range(c(mint,loV,upV))
    else
      xlim = range(c(mint,mV))

    if(is.null(ylim))
      if(!any(is.na(loV)))
        ylim = range(c(mint,loV,upV))
    else
      ylim = range(c(mint,mV))

    matplot(
      mint,
      mV,
      xlab = 'RMS(uE)',
      ylab = 'SD(E)',
      xlim = xlim,
      # xaxs = 'i',
      ylim = ylim,
      type = 'b',
      lty  = 3,
      pch  = 16,
      lwd  = lwd,
      cex  = ifelse(slide,0.5,1),
      col  = cols[col],
      main = title,
      log  = ifelse(logX,'xy','')
    )
    grid(equilogs = FALSE)
    abline(a = 0, b = 1,
           lty = 2,
           col = cols[1],
           lwd = lwd)

  } else {
    lines(
      mint, mV,
      type = 'b',
      lty  = 3,
      pch  = 16,
      lwd  = lwd,
      cex  = ifelse(slide,0.5,1),
      col  = cols[col]
    )

  }

  if(!any(is.na(loV)))
    if(slide) {
      ipl = seq(1,length(mint),length.out=nBin)
      polygon(
        c(mint,rev(mint)),
        c(loV, rev(upV)),
        col = cols_tr[col],
        border = NA)
      segments(
        mint[ipl], loV[ipl],
        mint[ipl], upV[ipl],
        col  = cols[col],
        lwd  = 1.5 * lwd,
        lend = 1)

    } else {
      segments(
        mint, loV,
        mint, upV,
        col  = cols[col],
        lwd  = 1.5 * lwd,
        lend = 1)

    }
  box()

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)


  invisible(
    list(
      mint   = mint,
      lwindx = intrv$lwindx,
      upindx = intrv$upindx,
      pc     = mV,
      pcl    = loV,
      pcu    = upV
    )
  )
}