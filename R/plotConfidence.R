#' Plot confidence curve for (uE,E) set
#'
#' @param E (vector) prediction uncertainty, uE, or predicted value, V
#' @param uE (vector) error or z-score
#' @param stat (function) statistic to use
#' @param oracle (logical) plot oracle curve
#' @param normal (logical) plot normal curve
#' @param conf_normal (logical) plot confidence band around normal
#' @param rep_normal (integer) sampling repetitions for normal curve
#' @param col (integer) color index for the curve
#' @param type (string) curve type: line ('l') or points ('p')
#' @param add (logical) add confidence curve to previous plot
#' @param xlab (string) x axis label
#' @param xlim (vector) limits of the x axis
#' @param ylab (string) y axis label
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param legend (string) legend for the dataset
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#'
#' @return Plot confidence curve for (E,uE) set

#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)  # Generate errors
#'   plotConfidence(E,uE)
#' }
#'
plotConfidence = function(
  E, uE,
  stat   = ErrViewLib::mue,
  oracle = TRUE,
  normal = FALSE,
  conf_normal = FALSE,
  rep_normal = 100,
  col    = 2,
  type   = c('l','p'),
  add    = FALSE,
  xlab   = 'k% discarded',
  xlim   = NULL,
  ylab   = 'MAE / MAE0',
  ylim   = NULL,
  title  = NULL,
  label  = 0,
  legend = NULL,
  gPars  = ErrViewLib::setgPars()
) {

  Sconf = function(X,pcVec=0:100) {
    S0 = stat(X)
    vstat = rep(0.0,length(pcVec))
    vstat[1] = 1.0
    for (i in 2:length(pcVec)) {
      k = pcVec[i]
      sel = 1:floor(k * M / 100)
      if(length(sel) == 0)
        vstat[i] = NA
      else
        vstat[i] = stat(X[-sel]) / S0
    }
    # Add 0% point
    return(vstat)
  }

  type = match.arg(type)

  if (as.numeric(length(E))*as.numeric(length(uE)) == 0)
    return()

  # Reorder data
  io  = order(uE, decreasing = TRUE)
  uE  = uE[io]
  E   = E[io]
  if(oracle) {
    io = order(abs(E), decreasing = TRUE) # Perfect set for oracle
    O  = E[io]
  }

  M = length(uE)
  pcVec = 0:100 # Vector of percentages

  vstat = Sconf(E,pcVec)

  if(oracle)
    vora = Sconf(O,pcVec)

  if(normal) {
    nrun = rep_normal
    vnorm = matrix(0, ncol = nrun, nrow = length(pcVec))
    for (i in 1:nrun) {
      X = rnorm(uE, 0, uE)
      vnorm[, i] = Sconf(X)
    }
    vnorm_mu = c()
    vnorm_mu[1] = 1
    if (conf_normal) {
      vnorm_lw = vnorm_up = c()
      vnorm_lw[1] = 1
      vnorm_up[1] = 1
    }
    for (i in 2:length(pcVec)) {
      X = vnorm[i, ]
      vnorm_mu[i] = mean(X)
      if (conf_normal) {
        ci = ErrViewLib::vhd(X)
        vnorm_lw[i] = ci[1]
        vnorm_up[i] = ci[2]
      }
    }
  }

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  if(add) {

    lines(pcVec, vstat,
          type = type,
          pch = 16,
          lwd = 2 * lwd,
          col = cols[col])

  } else {

    par(
      mfrow = c(1, 1),
      mar = mar,
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd,
      xaxs = 'i',
      yaxs = 'i',
      cex.main = 1
    )

    if (is.null(xlim))
      xlim = range(pcVec)

    if (is.null(ylim))
      ylim = c(0, 1)

    plot(
      pcVec, vstat,
      type = type,
      lty  = 1,
      pch  = 16,
      lwd  = 2*lwd,
      col  = cols[col],
      xlab = xlab,
      xlim = xlim,
      ylab = ylab,
      ylim = ylim,
      main = title
    )
    grid()

    if(oracle)
      lines(pcVec, vora, lty=2, lwd = 2*lwd, col=cols[1])

    if(normal) {
      lines(pcVec, vnorm_mu, lty = 3, lwd = 2*lwd, col=cols[1])
      if(conf_normal) {
        polygon(
          c(pcVec,rev(pcVec)),
          c(vnorm_up,rev(vnorm_lw)),
          col = cols_tr[1],
          border = NA)
      }
    }

    lty = if(type == 'l') 1 else NA
    pch = if(type == 'l') NA else 16
    lcol = col
    if(is.null(legend))
      legend = 'Data'
    if(oracle) {
      legend = c('Oracle',legend)
      lty = c(2,lty)
      pch = c(NA,pch)
      lcol = c(1,lcol)
    }
    if(normal) {
      legend = c('Normal',legend)
      lty = c(3,lty)
      pch = c(NA,pch)
      lcol = c(1,lcol)
    }
    legend(
      'bottomleft', bty = 'n', inset = 0.05,
      legend = legend,
      lty = lty,
      lwd = 2*lwd,
      col = cols[lcol],
      pch = pch
    )

    box()

    if(label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj = 1,
        cex = cex,
        line = 0.3)
  }
}