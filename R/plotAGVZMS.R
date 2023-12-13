#' Adversarial Group Validation of calibration based on the
#' fraction of validated confidence intervals of the <Z^2>
#' statistic
#'
#' @param Z (vector) set of z-score values to be tested
#' @param mZ2 (numeric) target value for <Z^2> (default `1`)
#' @param popMin (integer) minimal bin count in an interval
#' @param popMax (integer) maximal bin count in an interval
#' @param method (string) method used to estimate 95 percent CI on <Z^2>
#' @param BSmethod (string) bootstrap variant
#' @param nBoot (integer) number of bootstrap replicas
#' @param nMC (integer) number of random intervals per size
#' @param add (logical) add to previous graph ?
#' @param ylim (vector) limits of the y axis
#' @param col (integer) color index of main curve
#' @param verbose (logical) print progress messages
#' @param control (logical) estimate AGV for control sample (normal-standard)
#' @param colControl (integer) color index of control curve
#' @param gPars (list) graphical parameters
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#'
#' @return Invisibly returns a list of worst stats. Mainly used
#'   for its plotting side effect.
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   plotAGVZMS(E/uE, method = 'stud')
#' }
plotAGVZMS <- function(
  Z,
  popMin   = 100,
  popMax   = length(Z)/10,
  mZ2      = 1,
  nBoot    = 1500,
  method   = c('bootstrap','stud','auto'),
  BSmethod = c('bca','perc','basic'),
  nMC      = 1000,
  add      = FALSE,
  ylim     = NULL,
  col      = 6,
  verbose  = FALSE,
  control  = TRUE,
  colControl = 2,
  title    = '',
  label    = 0,
  gPars    = ErrViewLib::setgPars()
) {

  method   = match.arg(method)
  BSmethod = match.arg(BSmethod)

  # Generate control sample
  if(control)
    Zc = rnorm(Z)

  M = length(Z)
  sizes  = M / 2^{1:20}
  sizes  = round(sizes[sizes >= popMin & sizes <= popMax])

  pVal = plo = pup = c()
  pValc = ploc = pupc = c()
  minZS = rep(NA, length(sizes))
  maxZS = rep(NA, length(sizes))
  for (i in seq_along(sizes)) {
    OK = OKc = 0
    for (j in 1:nMC) {
      if(verbose & j%%10 ==0)
        cat(sizes[i],' : ',j,' / ',nMC,'\n')
      sam = sample(M, sizes[i])
      z = Z[sam]
      zs = ErrViewLib::ZMSCI(
        z, method = method,
        nBoot = nBoot, CImethod = BSmethod)
      ok = (zs$ci[1] - mZ2) * (zs$ci[2] - mZ2) <= 0
      OK = OK + as.numeric(ok)
      if(!ok) {
        minZS[i] = min(minZS[i],zs$mean, na.rm = TRUE)
        maxZS[i] = max(maxZS[i],zs$mean, na.rm = TRUE)
      }
      if(control) {
        z = Zc[sam]
        zs = ErrViewLib::ZMSCI(
          z, method = method,
          nBoot = nBoot, CImethod = BSmethod)
        ok = (zs$ci[1] - mZ2) * (zs$ci[2] - mZ2) <= 0
        OKc = OKc + as.numeric(ok)
      }
    }
    pVal[i] = OK / nMC
    ci      = DescTools::BinomCI(OK, nMC, method = "wilsoncc")
    plo[i]  = ci[,2]
    pup[i]  = ci[,3]
    if(control) {
      pValc[i] = OKc / nMC
      ci       = DescTools::BinomCI(OKc, nMC, method = "wilsoncc")
      ploc[i]  = ci[,2]
      pupc[i]  = ci[,3]
    }
  }

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

  if(!add) {
    if(is.null(ylim))
      ylim = c(min(plo),1)

    plot(
      sizes/M, pVal,
      type = 'b', pch = 19, cex = 0.75, col = cols[col],
      log = 'x', xlab = 'Relative group size',
      ylim = ylim, yaxs = 'i',
      ylab = 'Fraction of valid <Z^2> CIs',
      main = title
    )
    grid(equilogs = FALSE)
    abline(h = 0.95, col = cols[1], lty = 2)
    if(control)
      points(
        sizes/M, pValc,
        type = 'b', pch = 19, col = cols[colControl], cex = 0.75)
    box()
    if(label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj = 1,
        cex = cex,
        line = 0.3)
  } else {
    points(
      sizes/M, pVal,
      type = 'b', pch = 19, col = cols[col], cex = 0.75)
  }
  segments(sizes/M,plo,sizes/M,pup, col = cols[col], lwd = lwd)
  if(control)
    segments(sizes/M,ploc,sizes/M,pupc, col = cols[colControl], lwd = lwd)


  invisible(
    list(
      sizes = sizes,
      minZS = minZS,
      maxZS = maxZS
    )
  )

}
#' Adversarial Group Calibration based on the Scal statistic
#'
#' @param Z (vector) set of z-score values to be tested
#' @param popMin (integer) minimal bin count in an interval
#' @param method (string) method used to estimate 95 percent CI on <Z^2>
#' @param nGroup (integer) number random groups sampled to select worst one
#' @param nMC (integer) number of repeats for worst group selection
#' @param control (logical) estimate AGV for control sample (normal-standard)
#' @param colControl (integer) color index of control curve
#' @param dist (string) model error distribution to generate the control
#'    values. One of 'Normal' (default), 'Uniform', 'Normpn', Normp4', 'Laplace', 'Tn' or 'T4'
#' @param df (integer) degrees of freedom for distributions 'Normpn' and 'Tn'
#' @param add (logical) add to previous graph ?
#' @param ylim (vector) limits of the y axis
#' @param col (integer) color index of main curve
#' @param gPars (list) graphical parameters
#' @param title (string) a title to display above the plot
#' @param legend (logical) add a legend (default: TRUE) ?
#' @param label (integer) index of letter for subplot tag
#'
#' @return Nothing.
#'
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   plotAGV(E/uE)
#' }
plotAGV<- function(
  Z,
  popMin   = 50,
  nGroup   = 100,
  nMC      = 500,
  add      = FALSE,
  ylim     = NULL,
  col      = 6,
  control  = TRUE,
  colControl = 2,
  stat     = c('max','mean','median','q95'),
  dist     = c('Normal','Uniform','Normpn','Normp4','Laplace','T4','Student'),
  df       = 4,
  title    = '',
  label    = 0,
  legend   = TRUE,
  gPars    = ErrViewLib::setgPars()
) {

  dist = match.arg(dist)
  stat = match.arg(stat)

  if(control) {
    # Unit-variance distributions
    Normal = function(N, ...)
      rnorm(N)
    Student = function(N, df = df)
      rt(N, df = df) / sqrt(df/(df-2))
    T4 = function(N, df = 4)
      rt(N, df = df) / sqrt(df/(df-2))
    Uniform = function(N)
      runif(N, -sqrt(3), sqrt(3))
    Normpn = function(N, df = df)
      normalp::rnormp(N, p = df) / sqrt(df^(2/df)*gamma(3/df)/gamma(1/df))
    Laplace = function(N, df = 1)
      normalp::rnormp(N, p = df) / sqrt(df^(2/df)*gamma(3/df)/gamma(1/df))
    Normp4 = function(N, df = 4)
      normalp::rnormp(N, p = df) / sqrt(df^(2/df)*gamma(3/df)/gamma(1/df))
    # Get distribution function from arg name
    distFun = get(dist)
    # Generate control sample
    Zc = distFun(length(Z), df)
  }


  M = length(Z)
  sizes  = M / 2^{0:20}
  sizes  = round(sizes[sizes >= popMin & sizes <= M])

  ag = agu = agc = agcu =rep(NA, length(sizes))
  for (i in seq_along(sizes)) {
    wtz = wtzc = c()
    for (j in 1:nMC) {
      tz = tzc = c()
      for (k in 1:nGroup) {
        sam = sample(M, sizes[i])
        tz[k] = abs(log(mean(Z[sam]^2)))
        if(control)
          tzc[k] = abs(log(mean(Zc[sam]^2)))
      }
      wtz[j] = switch(
        stat,
        max    = max(tz),
        mean   = mean(tz),
        median = median(tz),
        q95    = quantile(tz, probs = 0.95)
      )
      if(control)
        wtzc[j] = switch(
          stat,
          max    = max(tzc),
          mean   = mean(tzc),
          median = median(tzc),
          q95    = quantile(tzc, probs = 0.95)
        )
    }
    ag[i]  = mean(wtz)
    agu[i] = sd(wtz)
    if(control) {
      agc[i]  = mean(wtzc)
      agcu[i] = sd(wtzc)
    }
  }

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
    cex.main = 1,
    xpd = FALSE # Clipping
  )
  x = sizes/M
  if(!add) {
    if(is.null(ylim))
      ylim = range(c(0,ag+2*agu,agc+2*agcu))
    plot(
      x, ag,
      type = 'b', pch = 19, cex = 0.75, col = cols[col],
      log = 'x', xlab = 'Relative group size',
      ylim = ylim, yaxs = 'i',
      ylab = 'Calibration error of worst group',
      main = title
    )
    grid(equilogs = FALSE)
    polygon(c(x,rev(x)),c(ag+2*agu,rev(ag-2*agu)),
            border = NA, col = cols_tr[col])
    if(control) {
      points(
        x, agc,
        type = 'b', pch = 19, col = cols[colControl], cex = 0.75)
      polygon(c(x,rev(x)),c(agc+2*agcu,rev(agc-2*agcu)),
              border = NA, col = cols_tr[colControl])
    }
    box()
    if(legend)
      legend(
        'topright', bty = 'n',
        legend = c('Data','Reference'),
        pch = 19, lty = 1, col = gPars$cols[c(6,2)])
    if(label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj = 1,
        cex = cex,
        line = 0.3)
  } else {
    points(
      x, ag,
      type = 'b', pch = 19, col = cols[col], cex = 0.75)
    polygon(c(x,rev(x)),c(ag+2*agu,rev(ag-2*agu)),
            border = NA, col = cols_tr[col])
  }

}