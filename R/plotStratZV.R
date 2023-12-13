#' Aggregate a stratified sample to avoid strata with counts smaller than a threshold
#'
#' @param X (vector) stratified sample to be aggregated
#' @param popMin (integer) minimal count in a stratum
#' @param greedy (logical) use greedy algorithm to merge strata (default: TRUE)
#'
#' @return a list containing the new X vector and the values and counts of the new strata
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Uncertainty
#'   X   = signif(uE,1)                # Stratify uE
#'   ST  = stratAgg(X)
#'   print(str(ST))
#' }
stratAgg <- function(X, popMin = 50, greedy = TRUE) {

  values = sort(unique(X)) #!!! Do not use names(table(X)) for this
                           # !!! (values are truncated)
  counts = as.vector(table(X))

  if (greedy) {
    # Merge strata until all counts > popMin
    while (sum(counts < popMin) != 0) {

      ## Stratum to merge
      first = which(counts < popMin)[1]

      # Select smallest neighbor target
      if(first == 1) {
        j = first + 1
      } else if(first == length(counts)) {
        j = first - 1
      } else {
        pm1 = counts[first - 1]
        pp1 = counts[first + 1]
        if(pm1 < pp1) {
          j = first - 1
        } else {
          j = first + 1
        }
      }

      ## Merge
      seli = X == values[first]
      selj = X == values[j]

      ### Weighted mean of merged strata
      values[first] =
        (counts[first] * values[first] +
           counts[j] * values[j]) /
        (counts[first] + counts[j])
      counts[first] = counts[first] + counts[j]
      counts[j] = 0
      X[seli] = values[first]
      X[selj] = values[first]

      sel    = counts != 0
      values = values[sel]
      counts = counts[sel]
    }

  } else {
    # Aggregate contiguous strata smaller than popMin
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

        # Weighted mean of merged strata
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
  }

  return(
    list(
      X      = X,
      values = values,
      counts = counts
    )
  )
}


#' Plot local z-score variance to assess calibration and tightness
#' for stratified conditioning variables
#'
#' @param X (vector) abscissae of the Z values
#' @param Z (vector) set of z-score values to be tested
#' @param aggregate (logical) aggregate contiguous strata smaller than popMin (default: TRUE)
#' @param popMin (integer) minimal count in a stratum
#' @param greedy (logical) use greedy algorithm to merge strata (default: TRUE)
#' @param varZ (numeric) target value for Var(Z) (default `1`)
#' @param method (string) method used to estimate 95 percent CI on Var(Z)
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
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
  greedy    = TRUE,
  plot      = TRUE,
  method    = c('cho', 'bootstrap', 'chisq','auto'),
  BSmethod  = c('bca', 'perc', 'basic'),
  nBoot     = 1500,
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
  if(aggregate) {
    ST     = stratAgg(X, popMin, greedy)
    X      = ST$X
    values = ST$values
    counts = ST$counts
  } else {
    values = sort(unique(X)) # !!! Do not use names(table(X)) for this
                             # !!! (values are truncated)
    counts = as.vector(table(X))
  }

  # Statistics
  mV = loV = upV  = c()
  for (i in seq_along(values)) {
    f   = values[i]
    sel = X == f
    c = ErrViewLib::varZCI(Z[sel], nBoot = nBoot,
                           method = method, CImethod = BSmethod)
    mV[i]  = c$mean
    loV[i] = c$ci[1]
    upV[i] = c$ci[2]
  }
  zs = varZCI(Z, method = 'auto', nBoot = nBoot, CImethod = BSmethod)
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
    colm = ifelse((loV0 - varZ) * (upV0 - varZ) < 0, cols[5], cols[2])
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
        line = 0.3
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
#' @param greedy (logical) use greedy algorithm to merge strata (default: TRUE)
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
  greedy    = TRUE,
  plot      = TRUE,
  xlab      = 'Conditioning variable',
  ylim      = NULL,
  title     = '',
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  N = length(Z)

  # Define strata
  if(aggregate) {
    ST     = stratAgg(X, popMin, greedy)
    X      = ST$X
    values = ST$values
    counts = ST$counts
  } else {
    values = sort(unique(X)) # !!! Do not use names(table(X)) for this
    # !!! (values are truncated)
    counts = as.vector(table(X))
  }

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
      ylab = '< Z >'
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
    # axis(3, at = xp[ss], labels = xp[ss])
    # mtext('Stratum index',
    #       side = 3,
    #       line = 2,
    #       cex = cex)

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
        line = 0.3
      )

  }

  # Statistics
  isVal   = loM * upM <= 0 & counts >= popMin
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
      pc     = mM,
      pcl    = loM,
      pcu    = upM,
      meanP  = mM0,
      meanPl = loM0,
      meanPu = upM0,
      fVal   = fVal,
      lofVal = lofVal,
      upfVal = upfVal,
      isVal  = isVal
    )
  )
}
#' Plot local <Z^2>
#' for stratified conditioning variables
#'
#' @param X (vector) conditioning variable
#' @param Z (vector) set of z-score values to be tested
#' @param mZ2 (numeric) target value for <Z^2> (default: 1)
#' @param method (string) method used to estimate 95 percent CI on <Z^2>
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
#' @param aggregate (logical) aggregate contiguous strata smaller than popMin
#' @param popMin (integer) minimal count in a stratum
#' @param greedy (logical) use greedy algorithm to merge strata (default: TRUE)
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
#'   plotStratZMS(X, E/uE, ylim = c(0, 2))
#' }
plotStratZMS = function(
  X, Z,
  mZ2       = 1,
  aggregate = TRUE,
  popMin    = 100,
  greedy    = TRUE,
  nBoot     = 1500,
  method    = c('bootstrap','stud','auto'),
  BSmethod  = c('bca','perc','basic'),
  plot      = TRUE,
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
  if(aggregate) {
    ST     = stratAgg(X, popMin, greedy)
    X      = ST$X
    values = ST$values
    counts = ST$counts
  } else {
    values = sort(unique(X)) # !!! Do not use names(table(X)) for this
    # !!! (values are truncated)
    counts = as.vector(table(X))
  }


  # LZM values
  mM = loM = upM = c()
  for (i in seq_along(values)) {
    sel    = X == values[i]
    zs     = ErrViewLib::ZMSCI(Z[sel], nBoot = nBoot,
                               method = method, CImethod = BSmethod )
    mM[i]   = zs$mean
    loM[i]  = zs$ci[1]
    upM[i]  = zs$ci[2]
  }
  zs    = ErrViewLib::ZMSCI(Z, method = 'auto', nBoot = nBoot,
                            CImethod = BSmethod)
  mM0   = zs$mean
  loM0  = zs$ci[1]
  upM0  = zs$ci[2]

  # Colors of symbols and segments

  if (length(gPars) == 0)
    gPars = ErrViewLib::setgPars()

  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  col = ifelse((loM - mZ2) * (upM-mZ2) <= 0, cols[5], cols[2])
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
      ylab = '< Z^2 >'
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
    # axis(3, at = xp[ss], labels = xp[ss])
    # mtext('Stratum index',
    #       side = 3,
    #       line = 2,
    #       cex = cex)

    abline(v = xp[ss], lwd = 1)
    segments(xp, loM, xp, upM, col = col, lwd = 2 * lwd)
    box()

    # Mean in margin
    ypos = par("usr")[4]
    pm   = round(mM0, digits = 2)
    colm = ifelse((loM0 - mZ2) * (upM0-mZ2) < 0, cols[5], cols[2])
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
    abline(h   = mZ2,
           lty = 2,
           col = cols[1],
           lwd = lwd)
    mtext(text = paste0(signif(mZ2,2),' -'),
          side = 2,
          at = mZ2,
          col = cols[1],
          cex = 0.75*cex,
          las = 1,
          adj = 1,
          font = 2)

    if (label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj = 1,
        cex = cex,
        line = 0.3
      )

  }

  # Statistics
  isVal   = (loM - mZ2) * (upM-mZ2) <= 0 & counts >= popMin
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
      pc     = mM,
      pcl    = loM,
      pcu    = upM,
      meanP  = mM0,
      meanPl = loM0,
      meanPu = upM0,
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
#' @param popMin (integer) minimal count in a stratum
#' @param greedy (logical) use greedy algorithm to merge strata (default: TRUE)
#' @param varZ (numeric) target value for Var(Z) (default `1`)
#' @param method (string) method used to estimate 95 percent CI on Var(Z)
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
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
  greedy    = TRUE,
  plot      = TRUE,
  nBoot     = 1500,
  xlab      = 'Conditioning variable',
  ylim      = NULL,
  title     = '',
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  # Define strata
  if(aggregate) {
    ST     = stratAgg(X, popMin, greedy)
    X      = ST$X
    values = ST$values
    counts = ST$counts
  } else {
    values = sort(unique(X)) # !!! Do not use names(table(X)) for this
    # !!! (values are truncated)
    counts = as.vector(table(X))
  }

  # Statistics
  mV = loV = upV  = c()
  for (i in seq_along(values)) {
    sel    = X == values[i]
    uEloc  = uE[sel]
    Eloc   = E[sel]
    bs     = boot::boot(cbind(uEloc, Eloc), ErrViewLib::rce, R = nBoot)
    bci    = boot::boot.ci(bs, conf = 0.95, type = 'bca')
    mV[i]  = bci$t0
    loV[i] = bci$bca[1, 4]
    upV[i] = bci$bca[1, 5]
  }
  bs   = boot::boot(cbind(uE, E), ErrViewLib::rce, R = nBoot)
  bci  = boot::boot.ci(bs, conf = 0.95, type = 'bca')
  mV0  = bci$t0
  loV0 = bci$bca[1, 4]
  upV0 = bci$bca[1, 5]

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
        line = 0.3
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