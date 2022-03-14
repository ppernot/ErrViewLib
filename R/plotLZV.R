#' Var(Var(Z)) by moments formula from Cho2008
#'
#' @param Z (vector) a data sample
#'
#' @return The variance on the sample variance
#' @export
#'
#' @examples
varvar = function (Z) {
  # Use formula (1) in CHo2008
  N = NROW(Z)
  mu = moments::all.moments(Z ,order.max=4, central=TRUE)[2:5]
  Vcho = (mu[4] - (N-3)/(N-1) * mu[2]^2 ) / N
  return(Vcho)
}
#' Estimate statistics of the variance of a sample
#'
#' @param Z (vector) a data sample
#' @param method (string) one of 'bootstrap' (default), 'cho' and
#'   'chisq' for the estimation of the statistics
#' @param CImethod (string) one of 'bca' (default), 'perc' and 'basic'
#'   for the CI estimation algorithm from bootstrap sample
#' @param nBoot (integer) number of bootstrap repeats
#' @param level (numeric) a probability level for CI (defauld: 0.95)
#' @param parallel (string) one of 'no' (default) and 'multicore'
#'
#' @return A list containing the mean, sd and ci for the variance
#'     of the Z sample
#' @export
#'
#' @examples
varZCI = function (
  Z,
  method = c('bootstrap','cho','chisq'),
  CImethod = c('bca','perc','basic'),
  nBoot = 1500,
  level = 0.95,
  parallel = c('no','multicore')
) {

  method   = match.arg(method)
  CImethod = match.arg(CImethod)
  parallel = match.arg(parallel)

  switch (
    method,
    bootstrap = {
      bst = boot::boot(Z, Hmisc::wtd.var, R = nBoot, stype = 'f',
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
      V  = var(Z)
      SD = sqrt(varvar(Z))
      list(
        mean   = V,
        sd     = SD,
        ci     = V + qnorm((1 + level) / 2) * c(-1, 1) * SD,
        method = method,
        level  = level
      )
    },
    chisq = {
      N = NROW(Z)
      V = var(Z)
      list(
        mean   = V,
        sd     = sqrt(2 / (N - 1)),
        ci     = c(
          qchisq((1 - level) / 2, df = N - 1) / (N - 1),
          qchisq((1 + level) / 2, df = N - 1) / (N - 1)
        ),
        method = method,
        level  = level
      )
    }
  )
}
#' Plot local Z-score variance to assess calibration and sharpness
#'
#' @param nBin  (integer) number of intervals for local coverage stats
#' @param ylim  (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot  (logical) plot the results
#' @param X     (vector) abscissae of the Z values
#' @param Z     (vector) set of Z values to be tested
#' @param logX  (logical) log-transform X
#' @param slide (logical) use sliding window for subsetting (X,Z)
#' @param method (string) method used to estimate 95% CI on Var(Z)
#' @param BSmethod (string) bootstrap variant, if method = 'bootstrap'
#' @param nBoot (integer) number of bootstrap replicas
#' @param xlab (string) X axis label
#' @param xlim (vector) min and max values of X axis
#'
#' @return Invisibly returns a list of LZV results. Mainly used
#'   for its plotting side effect.
#' @export
#'
#' @examples
plotLZV = function(
  X, Z,
  logX      = FALSE,
  nBin      = NULL,
  plot      = TRUE,
  slide     = NULL,
  method    = c('bootstrap','cho','chisq'),
  BSmethod  = c('bca','perc','basic'),
  nBoot     = 1500,
  xlab      = 'Calculated value',
  xlim      = NULL,
  ylim      = NULL,
  title     = '',
  label     = 0,
  gPars     = NULL
) {

  method   = match.arg(method)
  BSmethod = match.arg(BSmethod)

  N = length(Z)

  if(is.null(nBin))
    nBin  = max(min(floor(N/150),15),2)
  if(nBin <= 0)
    stop('>>> nBin should be > 0')
  if(is.null(slide))
    slide = nBin <= 4

  ord  = order(X)
  xOrd = X[ord]
  zOrd = Z[ord]

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

  # LZV values
  mV = loV = upV = mint = c()
  for (i in 1:nbr) {
    sel = lwindx[i]:upindx[i]
    M = length(sel)
    zLoc = zOrd[sel]
    zs = varZCI(zLoc, nBoot = nBoot, method = method, CImethod = BSmethod)
    mV[i]  = zs$mean
    loV[i] = zs$ci[1]
    upV[i] = zs$ci[2]
    mint[i] = 0.5*sum(range(xOrd[sel])) # Center of interval

  }
  zs = varZCI(zOrd, nBoot = nBoot, method = method, CImethod = BSmethod)
  mV0 = zs$mean
  loV0 = zs$ci[1]
  upV0 = zs$ci[2]

  if(plot) {
    # Plot ----
    if(length(gPars) == 0)
      gPars = ErrViewLib::setgPars()

    for (n in names(gPars))
      assign(n, rlist::list.extract(gPars, n))

    par(
      mfrow = c(1, 1),
      mar = c(mar[1:3],3),
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd
    )

    if(is.null(xlim))
      xlim = range(xOrd)

    if(is.null(ylim))
      ylim = range(c(loV, upV))

    matplot(
      mint,
      mV,
      xlab = xlab,
      ylab = 'Local Z-score variance',
      xlim = xlim,
      xaxs = 'i',
      ylim = ylim,
      type = 'b',
      lty = 3,
      pch = 16,
      lwd = lwd,
      cex = ifelse(slide,0.5,1),
      col  = cols[5],
      main = title,
      log = ifelse(logX,'x','')
    )
    grid()

    if(slide) {
      ipl = seq(1,length(mint),length.out=nBin)
      polygon(
        c(mint,rev(mint)),
        c(loV, rev(upV)),
        col = cols_tr[5],
        border = NA)
      segments(
        mint[ipl], loV[ipl],
        mint[ipl], upV[ipl],
        col  = cols[5],
        lwd  = 1.5 * lwd,
        lend = 1)

    } else {
      segments(
        mint, loV,
        mint, upV,
        col  = cols[5],
        lwd  = 1.5 * lwd,
        lend = 1)

    }
    xpos = pretty(xOrd)
    abline(h   = 1,
           lty = 2,
           col = cols[5],
           lwd = lwd)

    box()

    # Mean variance
    ypos = par("usr")[4]
    pm = signif(mV0,2)
    mtext(text = c(' Mean',paste0('- ',pm)),
          side = 4,
          at = c(ypos,mV0),
          col = c(1,cols[5]),
          cex = 0.75*cex,
          las = 1,
          font = 2)
    segments(
      xlim[2],loV0,
      xlim[2],upV0,
      col  = cols[5],
      lwd  = 6 * lwd,
      lend = 1
    )

    if(label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj = 1,
        cex = cex,
        line = 0.3)

  }

  invisible(
    list(
      mint   = mint,
      lwindx = lwindx,
      upindx = upindx,
      pc     = mV,
      pcl    = loV,
      pcu    = upV,
      meanP  = mV0
    )
  )
}
