#' Auxilliary function for plotConfidence
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
  pcVec = 0:99,
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
#' @param dfpr (logical) estimate Distance From Probabilistic Reference
#' @param col (integer) color index for the curve
#' @param type (string) curve type: line ('l') or points ('p')
#' @param add (logical) add confidence curve to previous plot
#' @param xlab (string) x axis label
#' @param xlim (vector) limits of the x axis
#' @param ylab (string) y axis label
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param showUk (logical) plot secondary axis with u_k values; supersedes `title`
#' @param showLegend (logical) display legend
#' @param legend (string) legend for the dataset
#' @param legLoc (string) location of legend (see \link[grDevices]{xy.coord})
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#'
#' @return Plot confidence curve for (E,uE) set
#'
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
  showUk = FALSE,
  dfpr   = TRUE,
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

  # Do we estimate area from probabilistic reference ?
  est_dfpr = dfpr & !normalize & probref

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
  pcVec = 0:99 # Vector of percentages

  # Statistics for data
  vstat = ErrViewLib::Sconf(E, pcVec, stat, normalize)

  # Statistics for oracle
  if(oracle)
    vora = ErrViewLib::Sconf(O, pcVec, stat, normalize)

  # Statistics for probabilistic reference
  # 1 - Draw samples of errors from uncertainties
  # 2 - Estimate Sconf over sample
  # 3 - Estimate mean and CI
  # 4 - Estimate dfpr
  vDFPR = vUP_DFPR  = NULL
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
    if(est_dfpr) {
      vDFPR = sum(abs(vstat-vnorm_mu), na.rm = TRUE)
      X = apply(vnorm, 2, function(x) sum(abs(x-vnorm_mu), na.rm = TRUE))
      vUP_DFPR = ErrViewLib::q95hd(X)
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
      main = ''
    )
    if(showUk) {
      axis(
        3,
        at = seq(20,80,by=20),
        labels = signif(
          quantile(uE,
                   probs = rev(seq(0.2,0.8,by = 0.2))),
          2),
        padj = 0.5
      )
      mtext(expression(u[k]),side = 3, at = 50, line = 1.5, cex = cex)
      # title(main = title, line = 2)

    } else {
      title(main = title)

    }
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
        lcol   = c(col,lcol)
      }
      if(est_dfpr) {
        legend = c(
          legend,
          paste0('DFPR = ',signif(vDFPR,2)),
          paste0('UP_95 = ',signif(vUP_DFPR,2))
        )
        lty = c(lty,NA,NA)
        pch = c(pch,NA,NA)
        lcol = c(lcol,1,1)
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

  invisible(
    list(
      DFPR  = vDFPR,
      UP_95 = vUP_DFPR
    )
  )
}
#' Plot confidence curve for (uE,E) set. Interface to plotConfidence
#' with different default arguments.
#'
#' @param E (vector) prediction uncertainty, uE, or predicted value, V
#' @param uE (vector) error or z-score
#' @param normalize (logical) use normalized statistic (default: FALSE)
#' @param statS (string) statistic to use. One of 'RMSE' (default), 'MAE' or 'Q95'
#' @param oracle (logical) plot oracle curve (default: FALSE)
#' @param probref (logical) plot probabilistic reference (probref) curve (default: TRUE)
#' @param conf_probref (logical) plot confidence band around probref curve (default: TRUE)
#' @param dist_probref (string) model error distribution to generate the reference curve. One of 'Normal' (default), 'Uniform', 'Normp4', 'Laplace' or 'T4'
#' @param rep_probref (integer) sampling repetitions for normal curve (default = 500)
#' @param dfpr (logical) estimate Distance From Probabilistic Reference
#' @param col (integer) color index for the curve (default: 6)
#' @param add (logical) add confidence curve to previous plot (default: FALSE)
#' @param xlab (string) x axis label (default: 'k \% discarded)
#' @param xlim (vector) limits of the x axis (default: NULL)
#' @param ylim (vector) limits of the y axis (default: NULL)
#' @param title (string) a title to display above the plot (default: '')
#' @param showUk (logical) plot secondary axis with u_k values; supersedes `title`
#' @param showLegend (logical) display legend (default: TRUE)
#' @param legend (string) legend for the dataset (default: NULL)
#' @param legLoc (string) location of legend (see \link[grDevices]{xy.coord}) (default: 'bottomleft')
#' @param label (integer) index of letter for subplot tag (default: 0)
#' @param gPars (list) graphical parameters (default: ErrViewLib::setgPars())
#'
#' @return Plot confidence curve for (E,uE) set
#'
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)  # Generate errors
#'   plotCC(E, uE, statS = 'MAE')
#' }
plotCC = function(
  E, uE,
  statS        = c('RMSE','MAE','Q95'),
  normalize    = FALSE,
  oracle       = FALSE,
  probref      = TRUE,
  conf_probref = TRUE,
  dist_probref = c('Normal','Uniform','Normp4','Laplace','T4'),
  rep_probref  = 500,
  dfpr         = TRUE,
  col          = 6,
  add          = FALSE,
  xlab         = 'k% discarded',
  xlim         = NULL,
  ylim         = NULL,
  title        = NULL,
  label        = 0,
  showUk       = FALSE,
  showLegend   = TRUE,
  legend       = NULL,
  legLoc       = 'bottomleft',
  gPars        = ErrViewLib::setgPars()
) {

  statS = match.arg(statS)
  dist_probref = match.arg(dist_probref)

  if(statS == 'RMSE') {
    stat = ErrViewLib::rmsd
    ylab = 'RMSE'
    if(normalize)
      ylab = 'RMSE / RMSE_0'
  } else if(statS == 'MAE') {
    stat = ErrViewLib::mue
    ylab = 'MAE'
    if(normalize)
      ylab = 'MAE / MAE_0'
  } else {
    stat = ErrViewLib::q95hd
    ylab = 'Q95'
    if(normalize)
      ylab = 'Q95 / Q95_0'
  }

  res = ErrViewLib::plotConfidence(
    E, uE,
    stat   = stat,
    normalize = normalize,
    oracle  = oracle,
    probref = probref,
    conf_probref = conf_probref,
    dist_probref = dist_probref,
    rep_probref = rep_probref,
    dfpr   = dfpr,
    col    = col,
    add    = add,
    xlab   = xlab,
    xlim   = xlim,
    ylab   = ylab,
    ylim   = ylim,
    title  = title,
    label  = label,
    showUk = showUk,
    showLegend = showLegend,
    legend = legend,
    legLoc = legLoc,
    gPars  = gPars
  )

  invisible(res)
}