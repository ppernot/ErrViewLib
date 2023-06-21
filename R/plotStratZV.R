#' Plot local z-score variance to assess calibration and tightness
#' for stratified conditioning variables
#'
#' @param X (vector) abscissae of the Z values
#' @param Z (vector) set of z-score values to be tested
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
#' @param score (logical) estimate calibration stats (default: `FALSE`)
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
  popMin    = 30,
  plot      = TRUE,
  method    = c('cho','bootstrap','chisq'),
  BSmethod  = c('bca','perc','basic'),
  nBoot     = 500,
  xlab      = 'Uncertainty',
  ylim      = NULL,
  title     = '',
  score     = FALSE,
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  method   = match.arg(method)
  BSmethod = match.arg(BSmethod)

  N = length(Z)

  # Define strata
  fH  = sort(unique(X))
  sel = table(X) > 1
  fH  = fH[sel]

  # Statistics
  mV = loV = upV = sH = c()
  for(i in seq_along(fH)) {
    f = fH[i]
    sel = X == f
    sH[i] = sum(sel)
    c = ErrViewLib::varZCI(Z[sel],
                           nBoot = max(nBoot, N),
                           method = method,
                           CImethod = BSmethod)
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

  if(length(gPars) == 0)
    gPars = ErrViewLib::setgPars()

  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  col = ifelse((loV-1)*(upV-1)<0, cols[5], cols[2])
  sel = sH < popMin
  col[sel] = 'gray'

  if(plot) {

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

    xp = 1:length(fH)
    xlim = c(0, length(xp)+1)

    if(is.null(ylim))
      ylim = range(c(loV, upV, varZ))

    plot(
      xp, mV,
      pch = 19, cex = 0.5, col = col,
      xaxt ='n', xlim = xlim, xaxs = 'i', xlab ='',
      ylim = ylim, yaxs = 'i', ylab = 'Var(Z)'
    )
    # grid()
    abline(h=1, lty = 2)

    # X axes
    ss = seq(1,max(xp),10)
    axis(1,at=xp[ss],labels = signif(fH[ss],1))
    mtext(xlab, side = 1, line = 2, cex = cex)
    axis(3,at=xp[ss],labels = xp[ss])
    mtext('Stratum index', side = 3, line = 2, cex = cex)

    abline(v=xp[ss],lwd=1)
    segments(xp,loV,xp,upV,col = col, lwd = 2*lwd)
    box()

    # Mean in margin
    ypos = par("usr")[4]
    pm   = round(mV0, digits = 2)
    colm = ifelse((loV0-1)*(upV0-1)<0, cols[5], cols[2])
    mtext(
      text = c(
        ' Mean',
        paste0('- ', pm)),
      side = 4,
      at   = c(ypos, mV0),
      col  = c(1,colm),
      cex  = 0.75*cex,
      las  = 1,
      font = 2)
    segments(
      xlim[2],loV0,
      xlim[2],upV0,
      col  = colm,
      lwd  = 6 * lwd,
      lend = 1
    )

    if(label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj = 1,
        cex = cex,
        line = 2)

  }

  fVal = NA
  if(score)
    fVal = sum(col == cols[5]) / sum(col != 'gray')

  invisible(
    list(
      fH     = fH,
      sH     = sH,
      pc     = mV,
      pcl    = loV,
      pcu    = upV,
      meanP  = mV0,
      meanPl = loV0,
      meanPu = upV0,
      fVal   = fVal
    )
  )
}
