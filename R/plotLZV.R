#' Var(Var(Z)) by moments formula from Cho2008
#'
#' @param Z (vector) a data sample
#'
#' @return The variance on the sample variance
#' @export
#'
varvar = function (Z) {
  # Use formula (1) in Cho2008
  N = NROW(Z)
  mu = moments::all.moments(Z ,order.max=4, central=TRUE)[2:5]
  return( (mu[4] - (N-3)/(N-1) * mu[2]^2 ) / N )
}
#' Estimate statistics of the variance of a sample
#'
#' @param Z (vector) a data sample
#' @param method (string) one of 'auto' (default), 'bootstrap', 'cho' and
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
  method   = c('bootstrap','cho','chisq','auto'),
  CImethod = c('bca','perc','basic'),
  nBoot = 1500,
  level = 0.95,
  parallel = c('no','multicore')
) {

  method   = match.arg(method)
  CImethod = match.arg(CImethod)
  parallel = match.arg(parallel)

  if(method == 'auto')
    method = ifelse(length(Z) >= 200, 'cho', 'bootstrap')

  switch (
    method,
    bootstrap = {
      bst = boot::boot(Z, Hmisc::wtd.var, R = max(nBoot,length(Z)), # for CIs
                       stype = 'f', parallel = parallel, normwt = TRUE)
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
#' Estimate statistics of <Z^2>
#'
#' @param Z (vector) a data sample
#' @param method (string) one of 'auto' (default), 'bootstrap' and 'stud'
#' @param CImethod (string) one of 'bca' (default), 'perc' and 'basic'
#'   for the CI estimation algorithm from bootstrap sample
#' @param nBoot (integer) number of bootstrap repeats
#' @param level (numeric) a probability level for CI (default: 0.95)
#' @param parallel (string) one of 'no' (default) and 'multicore'
#'
#' @return A list containing the mean, sd and ci for <Z^2>
#' @export
#'
ZMSCI = function (
  Z,
  method   = c('bootstrap','stud','auto','studc'),
  CImethod = c('bca','perc','basic'),
  nBoot = 1500,
  level = 0.95,
  parallel = c('no','multicore')
) {

  method   = match.arg(method)
  CImethod = match.arg(CImethod)
  parallel = match.arg(parallel)

  if(method == 'auto')
    method = ifelse(length(Z) >= 200, 'stud', 'bootstrap')

  switch (
    method,
    bootstrap = {
      bst = boot::boot(Z^2, ErrViewLib::mse,
                       R =  max(nBoot,length(Z)), # for CIs
                       parallel = parallel)
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
    stud = {
      mu = mean(Z^2)
      s2 = var(Z^2)
      N  = length(Z)
      sm = sqrt(s2/N)
      list(
        mean   = mu,
        sd     = sm,
        ci     = c(
          mu  + sm * qt((1 - level) / 2, df = N-1) ,
          mu  + sm * qt((1 + level) / 2, df = N-1)
        ),
        method = method,
        level  = level
      )
    },
    studc = {
      mu = mean(Z^2)
      s2 = var(Z^2)
      m3 = moments::moment(Z^2, order = 3, central = TRUE)
      N  = length(Z)
      sm = sqrt(s2/N)
      list(
        mean   = mu,
        sd     = sm,
        ci     = c(
          mu + m3/(6*s2*N) + sm * qt((1 - level) / 2, df = N-1) ,
          mu + m3/(6*s2*N) + sm * qt((1 + level) / 2, df = N-1)
        ),
        method = method,
        level  = level
      )
    }
  )
}
#' Plot local z-score variance to assess calibration and tightness
#'
#' @param X (vector) abscissae of the Z values
#' @param Z (vector) set of z-score values to be tested
#' @param aux (vector) auxilliary vector to resolve ties in X sorting
#' @param varZ (numeric) target value for Var(Z) (default `1`)
#' @param logX (logical) log-transform X
#' @param method (string) method used to estimate 95 percent CI on Var(Z)
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
#' @param intrv (object) intervals generated by `genIntervals` (default: `NULL`)
#' @param nBin (integer) number of intervals for local coverage stats
#' @param slide (logical) use sliding window for subsetting (X,Z)
#' @param equiPop (logical) generate intervals  with equal bin counts
#'   (default: `equiPop = TRUE`)
#' @param popMin (integer) minimal bin count in an interval
#' @param logBin (logical) if `equiPop = FALSE`, one can choose between
#'   equal range intervals, or equal log-range intervals
#'   (default `logBin = TRUE`)
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#' @param xlab (string) X axis label
#' @param xlim (vector) min and max values of X axis
#' @param score (logical) estimate calibration stats (default: `FALSE`)
#' @param add (logical) add to previous graph ?
#' @param col (integer) color index of curve to add
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
  aux       = NULL,
  varZ      = 1,
  logX      = FALSE,
  nBin      = NULL,
  equiPop   = TRUE,
  popMin    = 30,
  logBin    = TRUE,
  intrv     = NULL,
  plot      = TRUE,
  slide     = FALSE,
  method    = c('bootstrap','cho','chisq','auto'),
  BSmethod  = c('bca','perc','basic'),
  nBoot     = 1500,
  xlab      = 'Conditioning variable',
  xlim      = NULL,
  ylim      = NULL,
  title     = '',
  score     = FALSE,
  add       = FALSE,
  col       = 5,
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  method   = match.arg(method)
  BSmethod = match.arg(BSmethod)

  N = length(Z)

  if(is.null(aux))
    ord  = order(X)
  else
    ord  = order(X,aux)
  xOrd = X[ord]
  zOrd = Z[ord]

  # Design local areas
  if(is.null(intrv)) {
    # if(is.null(slide))
    #   slide = nBin <= 4

    if(is.null(nBin))
      nBin  = max(min(floor(N/150),15),2)
    if(nBin <= 0)
      stop('>>> nBin should be > 0')
    Xin = N
    if(!equiPop)
      Xin = xOrd
    intrv = ErrViewLib::genIntervals(Xin, nBin, slide, equiPop, popMin, logBin)
  }
  nBin0 = nBin # Used if slide = TRUE
  nBin  = intrv$nbr

  # LZV values
  mV = loV = upV = mint = c()
  for (i in 1:nBin) {
    sel  = intrv$lwindx[i]:intrv$upindx[i]
    zs   = ErrViewLib::varZCI( zOrd[sel], nBoot = nBoot,
                   method = method, CImethod = BSmethod)
    mV[i]   = zs$mean
    loV[i]  = zs$ci[1]
    upV[i]  = zs$ci[2]
    mint[i] = mean(range(xOrd[sel])) # Center of interval

  }
  zs = ErrViewLib::varZCI(zOrd, method = 'auto',
                          nBoot = nBoot, CImethod = BSmethod)
  mV0  = zs$mean
  loV0 = zs$ci[1]
  upV0 = zs$ci[2]

  colors = ifelse((loV-varZ)*(upV-varZ) <= 0, col, 2)

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
        ylim = range(c(loV, upV, varZ))

      plot(
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
        col  = cols[colors],
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
        col  = cols[colors]
      )
    }

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
        col  = cols[colors[ipl]],
        lwd  = 1.5 * lwd,
        lend = 1)

    } else {
      segments(
        mint, loV,
        mint, upV,
        col  = cols[colors],
        lwd  = 1.5 * lwd,
        lend = 1)

    }
    box()

    # Mean variance
    ypos = par("usr")[4]
    pm   = round(mV0, digits = 2)
    mtext(
      text = c(
        ' Average',
        paste0('- ', pm)),
      side = 4,
      at   = c(ypos, mV0),
      col  = c(1,cols[col]),
      cex  = 0.75*cex,
      las  = 1,
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

  ZVE = ZVEUp = ZVMs = ZVM =
    fVal = lofVal = upfVal = NA
  if(score) {
    scores = abs(log(mV))

    # Max deviation
    im = which.max(scores)
    ZVM = exp( sign(log(mV[im])) * scores[im] )

    # Significant ?
    ZVMs = FALSE
    if(varZ < loV[im] | varZ > upV[im])
      ZVMs = TRUE

    # Mean deviation
    ZVE = exp(mean(scores))
    scores = pmax(log(upV/mV),log(mV/loV))
    ZVEUp = exp(mean(scores))

    # Fraction of valid intervals
    success = sum((loV-varZ)*(upV-varZ) <= 0)
    trials  = length(loV)
    fVal    = success/trials
    ci      = DescTools::BinomCI(success, trials, method = "wilsoncc")
    lofVal  = ci[,2]
    upfVal  = ci[,3]
  }

  invisible(
    list(
      mint   = mint,
      lwindx = intrv$lwindx,
      upindx = intrv$upindx,
      pc     = mV,
      pcl    = loV,
      pcu    = upV,
      meanP  = mV0,
      meanPl = loV0,
      meanPu = upV0,
      ZVE    = ZVE,
      ZVEUp  = ZVEUp,
      ZVM    = ZVM,
      ZVMs   = ZVMs,
      fVal   = fVal,
      lofVal = lofVal,
      upfVal = upfVal,
      isVal  = (loV-varZ)*(upV-varZ) <= 0
    )
  )
}
#' Plot local z-score mean to assess unbiasedness
#'
#' @param X (vector) abscissae of the Z values
#' @param Z (vector) set of z-score values to be tested
#' @param aux (vector) auxilliary vector to resolve ties in X sorting
#' @param logX (logical) log-transform X
#' @param intrv (object) intervals generated by `genIntervals` (default: `NULL`)
#' @param nBin (integer) number of intervals for local coverage stats
#' @param slide (logical) use sliding window for subsetting (X,Z)
#' @param equiPop (logical) generate intervals  with equal bin counts
#'   (default: `equiPop = TRUE`)
#' @param popMin (integer) minimal bin count in an interval
#' @param logBin (logical) if `equiPop = FALSE`, one can choose between
#'   equal range intervals, or equal log-range intervals
#'   (default `logBin = TRUE`)
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#' @param xlab (string) X axis label
#' @param xlim (vector) min and max values of X axis
#' @param add (logical) add to previous graph ?
#' @param col (integer) color index of curve to add
#'
#' @return Invisibly returns a list of LZM results. Mainly used
#'   for its plotting side effect.
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   plotLZM(uE, E/uE, ylim = c(-1,1))
#' }
plotLZM = function(
  X, Z,
  aux       = NULL,
  logX      = FALSE,
  nBin      = NULL,
  equiPop   = TRUE,
  popMin    = 30,
  logBin    = TRUE,
  intrv     = NULL,
  plot      = TRUE,
  slide     = FALSE,
  xlab      = 'Conditioning variable',
  xlim      = NULL,
  ylim      = NULL,
  title     = '',
  add       = FALSE,
  col       = 5,
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  N = length(Z)

  if(is.null(aux))
    ord  = order(X)
  else
    ord  = order(X,aux)
  xOrd = X[ord]
  zOrd = Z[ord]

  # Design local areas
  if(is.null(intrv)) {
    # if(is.null(slide))
    #   slide = nBin <= 4

    if(is.null(nBin))
      nBin  = max(min(floor(N/150),15),2)
    if(nBin <= 0)
      stop('>>> nBin should be > 0')
    Xin = N
    if(!equiPop)
      Xin = xOrd
    intrv = ErrViewLib::genIntervals(Xin, nBin, slide, equiPop, popMin, logBin)
  }
  nBin0 = nBin # Used if slide = TRUE
  nBin  = intrv$nbr

  # LZM values
  mM = loM = upM = mint = c()
  for (i in 1:nBin) {
    sel    = intrv$lwindx[i]:intrv$upindx[i]
    M      = length(sel)
    zLoc   = zOrd[sel]
    mM[i]  = mean(zLoc)
    loM[i] = mM[i] + sd(zLoc)/sqrt(M) * qt(0.025, df = M-1)
    upM[i] = mM[i] + sd(zLoc)/sqrt(M) * qt(0.975, df = M-1)
    mint[i] = mean(range(xOrd[sel])) # Center of interval
  }
  mM0   = mean(Z)
  loM0  = mM0 - sd(Z)/sqrt(N) * 1.96
  upM0  = mM0 + sd(Z)/sqrt(N) * 1.96

  colors = ifelse(loM*upM <= 0, col, 2) # also used for scores...

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
        ylim = range(c(loM, upM, 0))

      plot(
        mint,
        mM,
        xlab = xlab,
        ylab = '< Z >',
        xlim = xlim,
        xaxs = 'i',
        ylim = ylim,
        type = 'b',
        lty  = 3,
        pch  = 16,
        lwd  = lwd,
        cex  = ifelse(slide,0.5,1),
        col  = cols[colors],
        main = title,
        log  = ifelse(logX,'x','')
      )
      grid(equilogs = FALSE)
      abline(h   = 0,
             lty = 2,
             col = cols[1],
             lwd = lwd)
      mtext(text = '0',
            side = 2,
            at   = 0,
            col  = cols[1],
            cex  = 0.75*cex,
            las  = 1,
            adj  = 1,
            font = 2)

    } else {
      lines(
        mint,
        mM,
        type = 'b',
        lty  = 3,
        pch  = 16,
        lwd  = lwd,
        cex  = ifelse(slide,0.5,1),
        col  = cols[colors]
      )
    }

    if(slide) {
      ipl = seq(1,length(mint),length.out=nBin0)
      polygon(
        c(mint,rev(mint)),
        c(loM, rev(upM)),
        col = cols_tr[col],
        border = NA)
      segments(
        mint[ipl], loM[ipl],
        mint[ipl], upM[ipl],
        col  = cols[colors[ipl]],
        lwd  = 1.5 * lwd,
        lend = 1)

    } else {
      segments(
        mint, loM,
        mint, upM,
        col  = cols[colors],
        lwd  = 1.5 * lwd,
        lend = 1)

    }
    box()

    # Mean variance
    ypos = par("usr")[4]
    pm   = round(mM0, digits = 2)
    mtext(
      text = c(
        ' Average',
        paste0('- ', pm)),
      side = 4,
      at   = c(ypos, mM0),
      col  = c(1,cols[col]),
      cex  = 0.75*cex,
      las  = 1,
      font = 2)
    segments(
      xlim[2],loM0,
      xlim[2],upM0,
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

  # Fraction of valid intervals
  success = sum(loM *upM <= 0)
  trials  = length(loM)
  fVal    = success / trials
  ci      = DescTools::BinomCI(success, trials, method = "wilsoncc")
  lofVal  = ci[,2]
  upfVal  = ci[,3]

  invisible(
    list(
      mint   = mint,
      lwindx = intrv$lwindx,
      upindx = intrv$upindx,
      pc     = mM,
      pcl    = loM,
      pcu    = upM,
      meanP  = mM0,
      meanPl = loM0,
      meanPu = upM0,
      fVal   = fVal,
      lofVal = lofVal,
      upfVal = upfVal,
      isVal  = loM * upM <= 0
    )
  )
}
#' Plot local RCE to assess calibration and tightness
#'
#' @param X (vector) abscissae of the Z values
#' @param Z (vector) set of z-score values to be tested
#' @param aux (vector) auxilliary vector to resolve ties in X sorting
#' @param logX (logical) log-transform X
#' @param nBoot (integer) number of bootstrap replicas
#' @param intrv (object) intervals generated by `genIntervals` (default: `NULL`)
#' @param nBin (integer) number of intervals for local coverage stats
#' @param slide (logical) use sliding window for subsetting (X,Z)
#' @param equiPop (logical) generate intervals  with equal bin counts
#'   (default: `equiPop = TRUE`)
#' @param popMin (integer) minimal bin count in an interval
#' @param logBin (logical) if `equiPop = FALSE`, one can choose between
#'   equal range intervals, or equal log-range intervals
#'   (default `logBin = TRUE`)
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#' @param xlab (string) X axis label
#' @param xlim (vector) min and max values of X axis
#' @param add (logical) add to previous graph ?
#' @param col (integer) color index of curve to add
#'
#' @return Invisibly returns a list of LZV results. Mainly used
#'   for its plotting side effect.
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   plotLRCE(uE, uE, E, ylim = c(-1, 1))
#' }
plotLRCE = function(
  X, uE, E,
  aux       = NULL,
  logX      = FALSE,
  nBin      = NULL,
  equiPop   = TRUE,
  popMin    = 30,
  logBin    = TRUE,
  intrv     = NULL,
  plot      = TRUE,
  slide     = FALSE,
  nBoot     = 500,
  xlab      = 'Conditioning variable',
  xlim      = NULL,
  ylim      = NULL,
  title     = '',
  add       = FALSE,
  col       = 5,
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  N = length(uE)

  if(is.null(aux))
    ord  = order(X)
  else
    ord  = order(X,aux)
  xOrd = X[ord]
  uOrd = uE[ord]
  eOrd = E[ord]

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
  nBin0 = nBin # Used if slide = TRUE
  nBin  = intrv$nbr

  # RCE values
  mV = loV = upV = mint = c()
  for (i in 1:nBin) {
    sel    = intrv$lwindx[i]:intrv$upindx[i]
    uEloc  = uOrd[sel]
    Eloc   = eOrd[sel]
    bs     = boot::boot(cbind(uEloc,Eloc), ErrViewLib::rce,
                        R=max(nBoot,length(sel)))
    bci    = boot::boot.ci(bs, conf = 0.95, type = 'bca' )
    mV[i]  = bci$t0
    loV[i] = bci$bca[1, 4]
    upV[i] = bci$bca[1, 5]
    mint[i] = mean(range(xOrd[sel])) # Center of interval
  }
  # Basic CI for full dataset (too long otherwise)
  bs   = boot::boot(cbind(uE, E), ErrViewLib::rce, R = nBoot)
  mV0  = bs$t0
  loV0 = quantile(bs$t,0.025)
  upV0 = quantile(bs$t,0.975)

  colors = ifelse(loV*upV <= 0, col, 2) # also used for scores...

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
        ylim = range(c(loV, upV, 0))

      plot(
        mint,
        mV,
        xlab = xlab,
        ylab = 'Relative Calibration Error',
        xlim = xlim,
        xaxs = 'i',
        ylim = ylim,
        type = 'b',
        lty  = 3,
        pch  = 16,
        lwd  = lwd,
        cex  = ifelse(slide,0.5,1),
        col  = cols[colors],
        main = title,
        log  = ifelse(logX,'x','')
      )
      grid(equilogs = FALSE)
      abline(h   = 0,
             lty = 2,
             col = cols[1],
             lwd = lwd)
      mtext(text = paste0(0,' -'),
            side = 2,
            at = 0,
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
        col  = cols[colors]
      )
    }

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
        col  = cols[colors[ipl]],
        lwd  = 1.5 * lwd,
        lend = 1)

    } else {
      segments(
        mint, loV,
        mint, upV,
        col  = cols[colors],
        lwd  = 1.5 * lwd,
        lend = 1)

    }
    box()

    # Mean variance
    ypos = par("usr")[4]
    pm   = round(mV0, digits = 2)
    mtext(
      text = c(
        ' Average',
        paste0('- ', pm)),
      side = 4,
      at   = c(ypos, mV0),
      col  = c(1,cols[col]),
      cex  = 0.75*cex,
      las  = 1,
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

  # Fraction of valid intervals
  fVal    = lofVal = upfVal = NA
  isVal   = loV * upV <= 0
  success = sum(isVal)
  trials  = length(loV)
  fVal    = success/trials
  ci      = DescTools::BinomCI(success, trials, method = "wilsoncc")
  lofVal  = ci[,2]
  upfVal  = ci[,3]

  invisible(
    list(
      mint   = mint,
      lwindx = intrv$lwindx,
      upindx = intrv$upindx,
      pc     = mV,
      pcl    = loV,
      pcu    = upV,
      meanP  = mV0,
      meanPl = loV0,
      meanPu = upV0,
      fVal   = fVal,
      lofVal = lofVal,
      upfVal = upfVal,
      isVal  = isVal
    )
  )
}
#' Plot local Mean Squared z-score <Z^2> to assess calibration and tightness
#'
#' @param X (vector) abscissae of the Z values
#' @param Z (vector) set of z-score values to be tested
#' @param aux (vector) auxilliary vector to resolve ties in X sorting
#' @param varZ (numeric) target value for Var(Z) (default `1`)
#' @param logX (logical) log-transform X
#' @param method (string) method used to estimate 95 percent CI on <Z^2>
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
#' @param intrv (object) intervals generated by `genIntervals` (default: `NULL`)
#' @param nBin (integer) number of intervals for local coverage stats
#' @param slide (logical) use sliding window for subsetting (X,Z)
#' @param equiPop (logical) generate intervals  with equal bin counts
#'   (default: `equiPop = TRUE`)
#' @param popMin (integer) minimal bin count in an interval
#' @param logBin (logical) if `equiPop = FALSE`, one can choose between
#'   equal range intervals, or equal log-range intervals
#'   (default `logBin = TRUE`)
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#' @param xlab (string) X axis label
#' @param xlim (vector) min and max values of X axis
#' @param score (logical) estimate calibration stats (default: `FALSE`)
#' @param add (logical) add to previous graph ?
#' @param col (integer) color index of curve to add
#'
#' @return Invisibly returns a list of LZMS results. Mainly used
#'   for its plotting side effect.
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   plotLZMS(uE, E/uE, method = 'cho', ylim = c(0,2))
#' }
plotLZMS = function(
  X, Z,
  aux       = NULL,
  varZ      = 1,
  logX      = FALSE,
  nBin      = NULL,
  equiPop   = TRUE,
  popMin    = 100,
  logBin    = TRUE,
  intrv     = NULL,
  plot      = TRUE,
  slide     = FALSE,
  nBoot     = 1500,
  method    = c('bootstrap','stud','auto'),
  BSmethod  = c('bca','perc','basic'),
  xlab      = 'Conditioning variable',
  xlim      = NULL,
  ylim      = NULL,
  title     = '',
  score     = FALSE,
  add       = FALSE,
  col       = 5,
  colInv    = 2,
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  method   = match.arg(method)
  BSmethod = match.arg(BSmethod)

  N = length(Z)

  if(is.null(aux))
    ord  = order(X)
  else
    ord  = order(X,aux)
  xOrd = X[ord]
  zOrd = Z[ord]

  # Design local areas
  if(is.null(intrv)) {
    # if(is.null(slide))
    #   slide = nBin <= 4

    if(is.null(nBin))
      nBin  = max(min(floor(N/150),15),2)
    if(nBin <= 0)
      stop('>>> nBin should be > 0')
    Xin = N
    if(!equiPop)
      Xin = xOrd
    intrv = ErrViewLib::genIntervals(Xin, nBin, slide, equiPop, popMin, logBin)
  }
  nBin0 = nBin # Used if slide = TRUE
  nBin  = intrv$nbr

  # LZV values
  mV = loV = upV = mint = c()
  for (i in 1:nBin) {
    sel  = intrv$lwindx[i]:intrv$upindx[i]
    zs   = ErrViewLib::ZMSCI(zOrd[sel], nBoot =nBoot,
                             method = method, CImethod = BSmethod )
    mV[i]   = zs$mean
    loV[i]  = zs$ci[1]
    upV[i]  = zs$ci[2]
    mint[i] = mean(range(xOrd[sel])) # Center of interval
  }
  zs = ErrViewLib::ZMSCI(Z, method = 'auto',
                         nBoot =nBoot, CImethod = BSmethod)
  mV0   = zs$mean
  loV0  = zs$ci[1]
  upV0  = zs$ci[2]

  colors = ifelse((loV-varZ)*(upV-varZ) <= 0, col, colInv)

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
        ylim = range(c(loV, upV, varZ))

      plot(
        mint,
        mV,
        xlab = xlab,
        ylab = '< Z^2 >',
        xlim = xlim,
        xaxs = 'i',
        ylim = ylim,
        type = 'b',
        lty  = 3,
        pch  = 16,
        lwd  = lwd,
        cex  = ifelse(slide,0.5,1),
        col  = cols[colors],
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
        col  = cols[colors]
      )
    }

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
        col  = cols[colors[ipl]],
        lwd  = 1.5 * lwd,
        lend = 1)

    } else {
      segments(
        mint, loV,
        mint, upV,
        col  = cols[colors],
        lwd  = 1.5 * lwd,
        lend = 1)

    }
    box()

    # Mean variance
    ypos = par("usr")[4]
    pm   = round(mV0, digits = 2)
    mtext(
      text = c(
        ' Average',
        paste0('- ', pm)),
      side = 4,
      at   = c(ypos, mV0),
      col  = c(1,cols[col]),
      cex  = 0.75*cex,
      las  = 1,
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

  ZVE = ZVEUp = ZVMs = ZVM =
    fVal = lofVal = upfVal = NA
  if(score) {
    scores = abs(log(mV))

    # Max deviation
    im = which.max(scores)
    ZVM = exp( sign(log(mV[im])) * scores[im] )

    # Significant ?
    ZVMs = FALSE
    if(varZ < loV[im] | varZ > upV[im])
      ZVMs = TRUE

    # Mean deviation
    ZVE = exp(mean(scores))
    scores = pmax(log(upV/mV),log(mV/loV))
    ZVEUp = exp(mean(scores))

    # Fraction of valid intervals
    success = sum((loV-varZ)*(upV-varZ) <= 0)
    trials  = length(loV)
    fVal    = success/trials
    ci      = DescTools::BinomCI(success, trials, method = "wilsoncc")
    lofVal  = ci[,2]
    upfVal  = ci[,3]
  }

  invisible(
    list(
      mint   = mint,
      lwindx = intrv$lwindx,
      upindx = intrv$upindx,
      pc     = mV,
      pcl    = loV,
      pcu    = upV,
      meanP  = mV0,
      meanPl = loV0,
      meanPu = upV0,
      ZVE    = ZVE,
      ZVEUp  = ZVEUp,
      ZVM    = ZVM,
      ZVMs   = ZVMs,
      fVal   = fVal,
      lofVal = lofVal,
      upfVal = upfVal,
      isVal  = (loV-varZ)*(upV-varZ) <= 0
    )
  )
}
