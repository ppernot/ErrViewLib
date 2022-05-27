#' Var(Var(Z)) by moments formula from Cho2008
#'
#' @param Z (vector) a data sample
#'
#' @return The variance on the sample variance
#' @export
#'
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
#' @param level (numeric) a probability level for CI (default: 0.95)
#' @param parallel (string) one of 'no' (default) and 'multicore'
#'
#' @return A list containing the mean, sd and ci for the variance
#'     of the Z sample
#' @export
#'
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
#' Plot local z-score variance to assess calibration and tightness
#'
#' @param nBin  (integer) number of intervals for local coverage stats
#' @param ylim  (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot  (logical) plot the results
#' @param X     (vector) abscissae of the Z values
#' @param Z     (vector) set of z-score values to be tested
#' @param logX  (logical) log-transform X
#' @param slide (logical) use sliding window for subsetting (X,Z)
#' @param method (string) method used to estimate 95 percent CI on Var(Z)
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
#' @param xlab (string) X axis label
#' @param xlim (vector) min and max values of X axis
#' @param add (logical) add to previous graph ?
#' @param col (interger) color index of curve to add
#'
#' @return Invisibly returns a list of LZV results. Mainly used
#'   for its plotting side effect.
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   plotLZV(uE, E/uE, method = 'cho', ylim = c(0,2))
#' }
plotLZV = function(
  X, Z,
  varZ      = 1,
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
  add       = FALSE,
  col       = 5,
  label     = 0,
  gPars     = ErrViewLib::setgPars()
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
  intrv = ErrViewLib::genIntervals(N, nBin, slide)

  # LZV values
  mV = loV = upV = mint = c()
  for (i in 1:intrv$nbr) {
    sel = intrv$lwindx[i]:intrv$upindx[i]
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
    if(length(gPars) == 0)
      gPars = ErrViewLib::setgPars()

    for (n in names(gPars))
      assign(n, rlist::list.extract(gPars, n))

    if(is.null(xlim))
      xlim = range(xOrd)

    if(!add) {
      par(
        mfrow = c(1, 1),
        mar = c(mar[1:3],3), # Reserve of right margin space
        mgp = mgp,
        pty = 's',
        tcl = tcl,
        cex = cex,
        lwd = lwd,
        cex.main = 1
      )


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
        lty  = 3,
        pch  = 16,
        lwd  = lwd,
        cex  = ifelse(slide,0.5,1),
        col  = cols[col],
        main = title,
        log  = ifelse(logX,'x','')
      )
      grid(equilogs = FALSE)
      abline(h   = varZ,
             lty = 2,
             col = cols[1],
             lwd = lwd)
      mtext(text = paste0(signif(varZ,2),' -'),
            side = 2,
            at = varZ,
            col = cols[1],
            cex = 0.75*cex,
            las = 1,
            adj = 1,
            font = 2)

    } else {
      lines(
        mint,
        mV,
        type = 'b',
        lty  = 3,
        pch  = 16,
        lwd  = lwd,
        cex  = ifelse(slide,0.5,1),
        col  = cols[col]
      )
    }

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
    xpos = pretty(xOrd)
    box()

    # Mean variance
    ypos = par("usr")[4]
    pm = signif(mV0,2)
    mtext(text = c(' Mean',paste0('- ',pm)),
          side = 4,
          at = c(ypos,mV0),
          col = c(1,cols[col]),
          cex = 0.75*cex,
          las = 1,
          font = 2)
    segments(
      xlim[2],loV0,
      xlim[2],upV0,
      col  = cols[col],
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
      lwindx = intrv$lwindx,
      upindx = intrv$upindx,
      pc     = mV,
      pcl    = loV,
      pcu    = upV,
      meanP  = mV0
    )
  )
}