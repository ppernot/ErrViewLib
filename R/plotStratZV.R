#' Plot local z-score variance to assess calibration and tightness
#' for stratified conditioning variables
#'
#' @param X (vector) abscissae of the Z values
#' @param Z (vector) set of z-score values to be tested
#' @param aggregate (logical) aggregate contiguous strata smaller than popMin
#' @param varZ (numeric) target value for Var(Z) (default `1`)
#' @param method (string) method used to estimate 95 percent CI on Var(Z)
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
#' @param popMin (integer) minimal count in a stratum
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#' @param xlab (string) X axis label
#'
#' @return Invisibly returns a list of LZV results. Mainly used
#'   for its plotting side effect.
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   X   = signif(uE,1)
#'   plotStratZV(X, E/uE, method = 'cho', ylim = c(0,2))
#' }
plotStratZV = function(
  X, Z,
  varZ      = 1,
  aggregate = TRUE,
  popMin    = 30,
  plot      = TRUE,
  method    = c('cho', 'bootstrap', 'chisq'),
  BSmethod  = c('bca', 'perc', 'basic'),
  nBoot     = 500,
  xlab      = 'Conditioning variable',
  ylim      = NULL,
  title     = '',
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  method   = match.arg(method)
  BSmethod = match.arg(BSmethod)

  N = length(Z)

  # Define strata
  values = sort(unique(X))
  counts = as.vector(table(X))

  if (aggregate)
    ## Aggregate contiguous strata smaller than popMin
    while (sum(counts < popMin) != 0) {
      sel = counts < popMin
      i = 1
      while (!sel[which(sel)[i] + 1]) {
        i = i + 1
        if (i > length(which(sel)))
          break
      }
      if (i > length(which(sel)))
        break
      first = which(sel)[i]
      j = first + 1
      while (sel[j]) {
        seli = X == values[first]
        selj = X == values[j]

        values[first] =
          (counts[first] * values[first] +
             counts[j] * values[j]) /
          (counts[first] + counts[j])
        counts[first] = counts[first] + counts[j]
        counts[j] = 0

        X[seli] = values[first]
        X[selj] = values[first]

        j = j + 1
        if (j > length(counts))
          break
      }
      sel = counts != 0
      values = values[sel]
      counts = counts[sel]
    }

  # Remove loners
  sel    = counts > 1
  values = values[sel]
  counts = counts[sel]

  # Statistics
  mV = loV = upV  = c()
  for (i in seq_along(values)) {
    f   = values[i]
    sel = X == f
    c = ErrViewLib::varZCI(
      Z[sel],
      nBoot = max(nBoot, sum(sel)),
      method = method,
      CImethod = BSmethod
    )
    mV[i]  = c$mean
    loV[i] = c$ci[1]
    upV[i] = c$ci[2]
  }
  zs = varZCI(Z,
              nBoot = max(nBoot, N),
              method = method,
              CImethod = BSmethod)
  mV0  = zs$mean
  loV0 = zs$ci[1]
  upV0 = zs$ci[2]

  # Colors of symbols and segments

  if (length(gPars) == 0)
    gPars = ErrViewLib::setgPars()

  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  col = ifelse((loV - varZ) * (upV - varZ) <= 0, cols[5], cols[2])
  sel = counts < popMin
  col[sel] = 'gray'

  if (plot) {
    par(
      mfrow = c(1, 1),
      mar = c(mar[1:3], 3),
      # Reserve of right margin space
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd,
      cex.main = 1
    )

    xp = 1:length(values)
    xlim = c(0, length(xp) + 1)

    if (is.null(ylim))
      ylim = range(c(loV, upV, varZ))

    plot(
      xp,
      mV,
      pch = 19,
      cex = 0.5,
      col = col,
      xaxt = 'n',
      xlim = xlim,
      xaxs = 'i',
      xlab = '',
      ylim = ylim,
      yaxs = 'i',
      ylab = 'Var(Z)'
    )
    # grid()
    abline(h = 1, lty = 2)

    # X axes
    ss = seq(1, max(xp), 10)
    axis(1, at = xp[ss], labels = prettyNum(values[ss], digits = 2))
    mtext(xlab,
          side = 1,
          line = 2,
          cex = cex)
    axis(3, at = xp[ss], labels = xp[ss])
    mtext('Stratum index',
          side = 3,
          line = 2,
          cex = cex)

    abline(v = xp[ss], lwd = 1)
    segments(xp, loV, xp, upV, col = col, lwd = 2 * lwd)
    box()

    # Mean in margin
    ypos = par("usr")[4]
    pm   = round(mV0, digits = 2)
    colm = ifelse((loV0 - 1) * (upV0 - 1) < 0, cols[5], cols[2])
    mtext(
      text = c(' Mean',
               paste0('- ', pm)),
      side = 4,
      at   = c(ypos, mV0),
      col  = c(1, colm),
      cex  = 0.75 * cex,
      las  = 1,
      font = 2
    )
    segments(
      xlim[2],
      loV0,
      xlim[2],
      upV0,
      col  = colm,
      lwd  = 6 * lwd,
      lend = 1
    )

    if (label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj = 1,
        cex = cex,
        line = 2
      )

  }

  # Statistics
  isVal   = (loV-varZ)*(upV-varZ) <= 0 & counts >= popMin
  success = sum(isVal)
  trials  = sum(counts >= popMin)
  fVal    = success / trials
  ci      = DescTools::BinomCI(success, trials, method = "wilsoncc")
  lofVal  = ci[, 2]
  upfVal  = ci[, 3]

  invisible(
    list(
      values = values,
      counts = counts,
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
#' Plot local z-score mean to assess unbiasedness
#' for stratified conditioning variables
#'
#' @param X (vector) abscissae of the Z values
#' @param Z (vector) set of z-score values to be tested
#' @param aggregate (logical) aggregate contiguous strata smaller than popMin
#' @param popMin (integer) minimal count in a stratum
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#' @param xlab (string) X axis label
#'
#' @return Invisibly returns a list of LZV results. Mainly used
#'   for its plotting side effect.
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   X   = signif(uE,1)
#'   plotStratZM(X, E/uE, ylim = c(-1,1))
#' }
plotStratZM = function(
  X, Z,
  aggregate = TRUE,
  popMin    = 30,
  plot      = TRUE,
  xlab      = 'Conditioning variable',
  ylim      = NULL,
  title     = '',
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  N = length(Z)

  # Define strata
  values = sort(unique(X))
  counts = as.vector(table(X))

  if (aggregate)
    ## Aggregate contiguous strata smaller than popMin
    while (sum(counts < popMin) != 0) {
      sel = counts < popMin
      i = 1
      while (!sel[which(sel)[i] + 1]) {
        i = i + 1
        if (i > length(which(sel)))
          break
      }
      if (i > length(which(sel)))
        break
      first = which(sel)[i]
      j = first + 1
      while (sel[j]) {
        seli = X == values[first]
        selj = X == values[j]

        values[first] =
          (counts[first] * values[first] +
             counts[j] * values[j]) /
          (counts[first] + counts[j])
        counts[first] = counts[first] + counts[j]
        counts[j] = 0

        X[seli] = values[first]
        X[selj] = values[first]

        j = j + 1
        if (j > length(counts))
          break
      }
      sel = counts != 0
      values = values[sel]
      counts = counts[sel]
    }

  # Remove loners
  sel    = counts > 1
  values = values[sel]
  counts = counts[sel]

  # LZM values
  mM = loM = upM = c()
  for (i in seq_along(values)) {
    sel    = X == values[i]
    M      = sum(sel)
    zLoc   = Z[sel]
    mZ     = mean(zLoc)
    sZ     = sd(zLoc) / sqrt(M)
    mM[i]  = mZ
    loM[i] = mZ + sZ * qt(0.025, df = M-1)
    upM[i] = mZ + sZ * qt(0.975, df = M-1)
  }
  mM0   = mean(Z)
  loM0  = mM0 - sd(Z)/sqrt(N) * 1.96
  upM0  = mM0 + sd(Z)/sqrt(N) * 1.96

  # Colors of symbols and segments

  if (length(gPars) == 0)
    gPars = ErrViewLib::setgPars()

  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  col = ifelse(loM * upM <= 0, cols[5], cols[2])
  sel = counts < popMin
  col[sel] = 'gray'

  if (plot) {
    par(
      mfrow = c(1, 1),
      mar = c(mar[1:3], 3),
      # Reserve of right margin space
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd,
      cex.main = 1
    )

    xp = 1:length(values)
    xlim = c(0, length(xp) + 1)

    if (is.null(ylim))
      ylim = range(c(loM, upM, 0))

    plot(
      xp,
      mM,
      pch = 19,
      cex = 0.5,
      col = col,
      xaxt = 'n',
      xlim = xlim,
      xaxs = 'i',
      xlab = '',
      ylim = ylim,
      yaxs = 'i',
      ylab = 'Mean(Z)'
    )
    # grid()
    abline(h = 0, lty = 2)

    # X axes
    ss = seq(1, max(xp), 10)
    axis(1, at = xp[ss], labels = prettyNum(values[ss], digits = 2))
    mtext(xlab,
          side = 1,
          line = 2,
          cex = cex)
    axis(3, at = xp[ss], labels = xp[ss])
    mtext('Stratum index',
          side = 3,
          line = 2,
          cex = cex)

    abline(v = xp[ss], lwd = 1)
    segments(xp, loM, xp, upM, col = col, lwd = 2 * lwd)
    box()

    # Mean in margin
    ypos = par("usr")[4]
    pm   = round(mM0, digits = 2)
    colm = ifelse(loM0 * upM0 < 0, cols[5], cols[2])
    mtext(
      text = c(' Mean',
               paste0('- ', pm)),
      side = 4,
      at   = c(ypos, mM0),
      col  = c(1, colm),
      cex  = 0.75 * cex,
      las  = 1,
      font = 2
    )
    segments(
      xlim[2],
      loM0,
      xlim[2],
      upM0,
      col  = colm,
      lwd  = 6 * lwd,
      lend = 1
    )

    if (label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj = 1,
        cex = cex,
        line = 2
      )

  }

  # Statistics
  isVal   = loV * upV <= 0 & counts >= popMin
  success = sum(isVal)
  trials  = sum(counts >= popMin)
  fVal    = success / trials
  ci      = DescTools::BinomCI(success, trials, method = "wilsoncc")
  lofVal  = ci[, 2]
  upfVal  = ci[, 3]

  invisible(
    list(
      values = values,
      counts = counts,
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
#' Plot local RCE to assess calibration and tightness
#' for stratified conditioning variables
#'
#' @param X (vector) abscissae of the Z values
#' @param Z (vector) set of z-score values to be tested
#' @param aggregate (logical) aggregate contiguous strata smaller than popMin
#' @param varZ (numeric) target value for Var(Z) (default `1`)
#' @param method (string) method used to estimate 95 percent CI on Var(Z)
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
#' @param popMin (integer) minimal count in a stratum
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#' @param xlab (string) X axis label
#'
#' @return Invisibly returns a list of RCE results. Mainly used
#'   for its plotting side effect.
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   X   = signif(uE,1)
#'   plotStratRCE(X, uE, E, ylim = c(-1,1))
#' }
plotStratRCE= function(
  X, uE, E,
  aggregate = TRUE,
  popMin    = 30,
  plot      = TRUE,
  nBoot     = 500,
  xlab      = 'Conditioning variable',
  ylim      = NULL,
  title     = '',
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  N = length(uE)

  # Define strata
  values = sort(unique(X))
  counts = as.vector(table(X))

  if (aggregate)
    ## Aggregate contiguous strata smaller than popMin
    while (sum(counts < popMin) != 0) {
      sel = counts < popMin
      i = 1
      while (!sel[which(sel)[i] + 1]) {
        i = i + 1
        if (i > length(which(sel)))
          break
      }
      if (i > length(which(sel)))
        break
      first = which(sel)[i]
      j = first + 1
      while (sel[j]) {
        seli = X == values[first]
        selj = X == values[j]

        values[first] =
          (counts[first] * values[first] +
             counts[j] * values[j]) /
          (counts[first] + counts[j])
        counts[first] = counts[first] + counts[j]
        counts[j] = 0

        X[seli] = values[first]
        X[selj] = values[first]

        j = j + 1
        if (j > length(counts))
          break
      }
      sel = counts != 0
      values = values[sel]
      counts = counts[sel]
    }

  # Remove loners
  sel    = counts > 1
  values = values[sel]
  counts = counts[sel]

  # Statistics
  fRCE = function(X, indx = 1:nrow(X)) {
    RMV    = sqrt(mean(X[indx,1]^2))
    RMSE   = sqrt(mean(X[indx,2]^2))
    return( (RMV - RMSE) / RMV )
  }

  mV = loV = upV  = c()
  for (i in seq_along(values)) {
    sel    = X == values[i]
    M      = sum(sel)
    uEloc  = uE[sel]
    Eloc   = E[sel]
    bs     = boot::boot(cbind(uEloc,Eloc), fRCE, R=nBoot)
    R      = bs$t0
    uR     = sd(bs$t)
    mV[i]  = R
    loV[i] = R + uR * qt(0.025, df = M - 1)
    upV[i] = R + uR * qt(0.975, df = M - 1)
  }
  bs   = boot::boot(cbind(uE, E), fRCE, R=nBoot)
  R    = bs$t0
  uR   = sd(bs$t)
  mV0  = R
  loV0 = R - uR * 1.96
  upV0 = R + uR * 1.96

  # Colors of symbols and segments

  if (length(gPars) == 0)
    gPars = ErrViewLib::setgPars()

  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  col = ifelse(loV * upV <= 0, cols[5], cols[2])
  sel = counts < popMin
  col[sel] = 'gray'

  if (plot) {
    par(
      mfrow = c(1, 1),
      mar = c(mar[1:3], 3),
      # Reserve of right margin space
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd,
      cex.main = 1
    )

    xp = 1:length(values)
    xlim = c(0, length(xp) + 1)

    if (is.null(ylim))
      ylim = range(c(loV, upV, 0))

    plot(
      xp,
      mV,
      pch = 19,
      cex = 0.5,
      col = col,
      xaxt = 'n',
      xlim = xlim,
      xaxs = 'i',
      xlab = '',
      ylim = ylim,
      yaxs = 'i',
      ylab = 'RCE = (RMV - RMSE) / RMV'
    )
    # grid()
    abline(h = 0, lty = 2)

    # X axes
    ss = seq(1, max(xp), 10)
    axis(1, at = xp[ss], labels = prettyNum(values[ss], digits = 2))
    mtext(xlab,
          side = 1,
          line = 2,
          cex = cex)
    axis(3, at = xp[ss], labels = xp[ss])
    mtext('Stratum index',
          side = 3,
          line = 2,
          cex = cex)

    abline(v = xp[ss], lwd = 1)
    segments(xp, loV, xp, upV, col = col, lwd = 2 * lwd)
    box()

    # Mean in margin
    ypos = par("usr")[4]
    pm   = round(mV0, digits = 2)
    colm = ifelse(loV0 * upV0 < 0, cols[5], cols[2])
    mtext(
      text = c(' Mean',
               paste0('- ', pm)),
      side = 4,
      at   = c(ypos, mV0),
      col  = c(1, colm),
      cex  = 0.75 * cex,
      las  = 1,
      font = 2
    )
    segments(
      xlim[2],
      loV0,
      xlim[2],
      upV0,
      col  = colm,
      lwd  = 6 * lwd,
      lend = 1
    )

    if (label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj = 1,
        cex = cex,
        line = 2
      )

  }

  # Statistics
  isVal   = loV * upV <= 0 & counts >= popMin
  success = sum(isVal)
  trials  = sum(counts >= popMin)
  fVal    = success / trials
  ci      = DescTools::BinomCI(success, trials, method = "wilsoncc")
  lofVal  = ci[, 2]
  upfVal  = ci[, 3]

  invisible(
    list(
      values = values,
      counts = counts,
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