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
      if(verbose)
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