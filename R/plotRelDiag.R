#' Estimate statistics of the SD of a sample
#'
#' @param Z (vector) a data sample
#' @param method (string) one of 'bootstrap' (default) and 'cho'
#' @param CImethod (string) one of 'bca' (default), 'perc' and 'basic'
#'   for the CI estimation algorithm from bootstrap sample
#' @param nBoot (integer) number of bootstrap repeats
#' @param level (numeric) a probability level for CI (default: 0.95)
#' @param parallel (string) one of 'no' (default) and 'multicore'
#'
#' @return A list containing the mean, sd and ci for the SD
#'     of the Z sample
#' @export
#'
sdZCI = function (
  Z,
  method = c('bootstrap','cho'),
  CImethod = c('bca','perc','basic'),
  nBoot = 1500,
  level = 0.95,
  parallel = c('no','multicore')
) {

  method   = match.arg(method)
  CImethod = match.arg(CImethod)
  parallel = match.arg(parallel)

  wtd.sd = function(x,weights=NULL,normwt=TRUE)
    sqrt(Hmisc::wtd.var(x,weights=weights,normwt=normwt))

  switch (
    method,
    bootstrap = {
      bst = boot::boot(Z, wtd.sd, R = nBoot, stype = 'f',
                       parallel = parallel, normwt = TRUE)
      bci = boot::boot.ci(bst, conf = level, type = CImethod )
      ci = switch(
        CImethod,
        bca   = bci$bca[1, 4:5],
        perc  = bci$perc[1, 4:5],
        basic = bci$basic[1, 4:5]
      )
      list(
        mean   = bci$t0,
        sd     = sd(bst$t),
        ci     = ci,
        method = method,
        CImethod = CImethod,
        level  = level
      )
    },
    cho = {
      # LUP formula for y = sqrt(Var(X))
      V  = var(Z)
      S  = sqrt(V)
      SD = sqrt(ErrViewLib::varvar(Z)) / (2*V^(1/2))
      list(
        mean   = S,
        sd     = SD,
        ci     = S + qnorm((1 + level) / 2) * c(-1, 1) * SD,
        method = method,
        level  = level
      )
    }
  )
}
#' Plot reliability diagram
#'
#' @param X     (vector) vector of uncertainties
#' @param Y     (vector) vector of errors
#' @param nBin  (integer) number of intervals for local coverage stats
#' @param slide (logical) use sliding window for subsetting (X,Z)
#' @param equiPop (logical) generate intervals  with equal bin counts
#'   (default: `equiPop = TRUE`)
#' @param popMin (integer) minimal bin count in an interval
#' @param logBin (logical) if `equiPop = FALSE`, one can choose between
#'   equal range intervals, or equal log-range intervals
#'   (default `logBin = TRUE`)
#' @param method (string) one of 'bootstrap' (default) and 'cho'
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
#' @param xlim (vector) min and max values of X axis
#' @param logX  (logical) log-transform X
#' @param ylim  (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param add (logical) add to previous graph ?
#' @param col (integer) color index of curve to add
#'
#' @return Used for its plotting side effect.
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)  # Generate errors
#'   plotRelDiag(uE, E, nBin = 6, nBoot = 500)
#' }
plotRelDiag = function(
  X, Y,
  logX      = FALSE,
  nBin      = NULL,
  equiPop   = TRUE,
  popMin    = 30,
  logBin    = TRUE,
  intrv     = NULL,
  method    = c('bootstrap','cho'),
  BSmethod  = c('bca','perc','basic'),
  nBoot     = 500,
  slide     = FALSE,
  xlim      = NULL,
  ylim      = NULL,
  title     = '',
  label     = 0,
  add       = FALSE,
  col       = 5,
  gPars     = ErrViewLib::setgPars()
) {

  method   = match.arg(method)
  BSmethod = match.arg(BSmethod)

  N = length(Y)

  ord  = order(X)
  xOrd = X[ord]
  yOrd = Y[ord]

  # Design local areas
  if(is.null(intrv)) {
    if(is.null(nBin))
      nBin  = max(min(floor(N/150),15),2)
    if(nBin <= 0)
      stop('>>> nBin should be > 0')
    Xin = N
    if(!equiPop)
      Xin = xOrd
    intrv = ErrViewLib::genIntervals(Xin, nBin, slide, equiPop, popMin, logBin)
  }
  nBin0 = nBin
  nBin = intrv$nbr

  mV = loV = upV = mint = c()
  for (i in 1:nBin) {
    sel  = intrv$lwindx[i]:intrv$upindx[i]
    M    = length(sel)
    zLoc = yOrd[sel]
    zs   = sdZCI(
      zLoc,
      nBoot = max(nBoot, M),
      method = method,
      CImethod = BSmethod
    )
    mV[i]   = zs$mean
    loV[i]  = zs$ci[1]
    upV[i]  = zs$ci[2]
    mint[i] = sqrt(mean(xOrd[sel] ^ 2))
  }

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
      ipl = seq(1,length(mint),length.out=nBin0)
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