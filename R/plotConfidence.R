#' Auxilliary function for plotConfidence²
#'
#' @param X (vector) errors
#' @param pcVec (vector) percentages
#' @param stat (function) statistiv
#' @param normalize (logical) normalize statistic
#'
#' @return Vector of CC stats
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)  # Generate errors
#'   plotConfidence(E,uE)
#' }
Sconf = function(
  X,
  pcVec = 0:100,
  stat = ErrViewLib::mue,
  normalize = TRUE
) {
  M = length(X)
  S0 = stat(X)
  vstat = rep(0.0, length(pcVec))
  vstat[1] = S0
  if(normalize)
    vstat[1] = 1
  for (i in 2:length(pcVec)) {
    k = pcVec[i]
    sel = 1:floor(k * M / 100)
    if (length(sel) == 0) {
      vstat[i] = NA
    } else {
      vstat[i] = stat(X[-sel])
      if(normalize)
        vstat[i] = vstat[i] / S0
    }

  }
  return(vstat)
}
#' Plot confidence curve for (uE,E) set
#'
#' @param E (vector) prediction uncertainty, uE, or predicted value, V
#' @param uE (vector) error or z-score
#' @param normalize (logical) use normalized statistic
#' @param stat (function) statistic to use
#' @param oracle (logical) plot Oracle curve
#' @param probref (logical) plot probabilistic reference (probref) curve
#' @param conf_probref (logical) plot confidence band around probref curve
#' @param dist_probref (string) model error distribution to generate Ideal curve
#' @param rep_probref (integer) sampling repetitions for normal curve
#' @param col (integer) color index for the curve
#' @param type (string) curve type: line ('l') or points ('p')
#' @param add (logical) add confidence curve to previous plot
#' @param xlab (string) x axis label
#' @param xlim (vector) limits of the x axis
#' @param ylab (string) y axis label
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param showLegend (logical) display legend
#' @param legend (string) legend for the dataset
#' @param legLoc (string) location of legend (see \link[grDevices]{xy.coord})
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
  normalize = TRUE,
  oracle  = TRUE,
  probref = FALSE,
  conf_probref = FALSE,
  dist_probref = 'Normal',
  rep_probref = 100,
  col    = 2,
  type   = c('l','p'),
  add    = FALSE,
  xlab   = 'k% discarded',
  xlim   = NULL,
  ylab   = ifelse(normalize,'MAE / MAE0','MAE'),
  ylim   = NULL,
  title  = NULL,
  label  = 0,
  showLegend = TRUE,
  legend = NULL,
  legLoc = 'bottomleft',
  gPars  = ErrViewLib::setgPars()
) {

  type = match.arg(type)

  if (as.numeric(length(E))*as.numeric(length(uE)) == 0)
    return()

  # Unit-variance distributions
  Normal = function(N)
    rnorm(N)
  T4 = function(N,df=4)
    rt(N, df = df) / sqrt(df/(df-2))
  Uniform = function(N)
    runif(N, -sqrt(3), sqrt(3))
  Laplace = function(N, df = 1)
    normalp::rnormp(N, p = df) / sqrt(df^(2/df)*gamma(3/df)/gamma(1/df))
  Normp4 = function(N, df = 4)
    normalp::rnormp(N, p = df) / sqrt(df^(2/df)*gamma(3/df)/gamma(1/df))

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

  # Statistics for data
  vstat = ErrViewLib::Sconf(E, pcVec, stat, normalize)

  # Statistics for oracle
  if(oracle)
    vora = ErrViewLib::Sconf(O, pcVec, stat, normalize)

  # Statistics for probabilistic reference
  # 1 - Draw samples of errors from uncertainties
  # 2 - Estimate Sconf over sample
  # 3 - Estimate mean and CI
  if(probref) {
    nrun = rep_probref
    vnorm = matrix(0, ncol = nrun, nrow = length(pcVec))
    fun = get(dist_probref)
    for (i in 1:nrun) {
      X = uE * fun(M)
      vnorm[, i] = ErrViewLib::Sconf(X, pcVec, stat, normalize)
    }
    vnorm_mu = apply(vnorm, 1, mean, na.rm = TRUE)
    if (conf_probref) {
      ci = apply(vnorm, 1, ErrViewLib::vhd)
      vnorm_lw = ci[1,]
      vnorm_up = ci[2,]
    }
  }

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  if(add) {

    if(probref) {
      lines(pcVec, vnorm_mu, lty = 3, lwd = 2*lwd, col=cols[col])
      if(conf_probref) {
        polygon(
          c(pcVec,rev(pcVec)),
          c(vnorm_up,rev(vnorm_lw)),
          col = cols_tr[col],
          border = NA)
      }
    }
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

    if (is.null(ylim)) {
      vs = vstat[is.finite(vstat)]
      ylim = c(0, max(vs))
      if(oracle)
        ylim = range(c(ylim,vora[is.finite(vora)]))
      if(probref)
        ylim = range(c(ylim,vnorm_mu[is.finite(vnorm_mu)]))
    }

    plot(
      pcVec, vstat,
      type = 'n',
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

    if(normalize)
      abline(h = 1, lwd = 2 * lwd, col = cols[1])

    if(oracle)
      lines(pcVec, vora, lty = 2, lwd = 2*lwd, col=cols[1])

    if(probref) {
      lines(pcVec, vnorm_mu, lty = 3, lwd = 2*lwd, col=cols[col])
      if(conf_probref) {
        polygon(
          c(pcVec,rev(pcVec)),
          c(vnorm_up,rev(vnorm_lw)),
          col = cols_tr[col],
          border = NA
        )
      }
    }

    lines(
      pcVec, vstat,
      type = type,
      lty  = 1,
      pch  = 16,
      lwd  = 2*lwd,
      col  = cols[col]
    )

    # Build and display legend
    if(showLegend) {
      lty = if(type == 'l') 1 else NA
      pch = if(type == 'l') NA else 16
      lcol = col

      if(is.null(legend))
        legend = 'Data'

      if(oracle) {
        legend = c('Oracle',legend)
        lty    = c(2,lty)
        pch    = c(NA,pch)
        lcol   = c(1,lcol)
      }
      if(probref) {
        legend = c('Prob. ref.',legend)
        lty    = c(3,lty)
        pch    = c(NA,pch)
        lcol   = c(1,lcol)
      }
      legend(
        legLoc, bty = 'n', inset = 0.05,
        legend = legend,
        lty    = lty,
        lwd    = 2*lwd,
        col    = cols[lcol],
        pch    = pch
      )
    }

    box()

    # Display subfigure label
    if(label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj = 1,
        cex = cex,
        line = 0.3)
  }
}