#' Auxilliary function for plotConfidence
#'
#' @param X (vector) errors
#' @param pcVec (vector) percentages
#' @param stat (function) statistic
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
#' @param plot (function) produce plot ?
#' @param confStat (logical) estimate DFPR and/or AUCO ?
#' @param score (logical) estimate DFPR and/or AUCO ?
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
  plot   = TRUE,
  score  = TRUE,
  confStat = score,
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
  est_dfpr = confStat & !normalize & probref
  est_auco = confStat & !normalize & oracle

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
  io = order(uE, decreasing = TRUE)
  uE = uE[io]
  E  = E[io]
  if(oracle) {
    # Perfect ordered set for oracle
    io = order(abs(E), decreasing = TRUE)
    O  = E[io]
  }

  M = length(uE)
  pcVec = 0:99 # Vector of percentages

  # Statistics for data
  vstat = ErrViewLib::Sconf(E, pcVec, stat, normalize)

  # Statistics for oracle
  vAUCO = NULL
  if(oracle) {
    vora = ErrViewLib::Sconf(O, pcVec, stat, normalize)
    if(est_auco)
      vAUCO = sum(abs(vstat-vora), na.rm = TRUE)
  }

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

  if(plot) {
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
        if(est_auco) {
          legend = c(
            legend,
            paste0('AUCO = ',signif(vAUCO,2))
          )
          lty = c(lty,NA)
          pch = c(pch,NA)
          lcol = c(lcol,1)
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
            paste0('UP95 = ',signif(vUP_DFPR,2))
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

  }

  invisible(
    list(
      DFPR = vDFPR,
      UP95 = vUP_DFPR,
      AUCO = vAUCO
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
#' @param plot (function) produce plot ?
#' @param confStat (logical) estimate DFPR and/or AUCO
#' @param score (logical) estimate DFPR and/or AUCO
#' @param col (integer) color index for the curve (default: 6)
#' @param add (logical) add confidence curve to previous plot (default: FALSE)
#' @param xlab (string) x axis label (default: 'k \% discarded)
#' @param xlim (vector) limits of the x axis (default: NULL)
#' @param ylim (vector) limits of the y axis (default: NULL)
#' @param title (string) a title to display above the plot (default: '')
#' @param showUk (logical) plot secondary axis with u_k values; supersedes `title`
#' @param showLegend (logical) display legend (default: TRUE)
#' @param legend (string) legend for the dataset (default: NULL)
#' @param unit (string) unit string to add to ylab
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
  statS        = c('RMSE','RMSD','MAE','Q95'),
  normalize    = FALSE,
  oracle       = FALSE,
  probref      = TRUE,
  conf_probref = TRUE,
  dist_probref = c('Normal','Uniform','Normp4','Laplace','T4'),
  rep_probref  = 500,
  plot         = TRUE,
  score        = TRUE,
  confStat     = score,
  col          = 6,
  add          = FALSE,
  xlab         = 'k% discarded',
  xlim         = NULL,
  ylim         = NULL,
  title        = NULL,
  label        = 0,
  showUk       = FALSE,
  showLegend   = TRUE,
  unit         = '',
  legend       = NULL,
  legLoc       = 'bottomleft',
  gPars        = ErrViewLib::setgPars()
) {

  statS = match.arg(statS)
  dist_probref = match.arg(dist_probref)

  if(statS == 'RMSE') {
    stat = ErrViewLib::rmse
    ylab = paste0('RMSE ',unit)
    if(normalize)
      ylab = 'RMSE / RMSE_0'
  } else if(statS == 'RMSD') {
    stat = ErrViewLib::rmsd
    ylab = paste0('RMSD ',unit)
    if(normalize)
      ylab = 'RMSD / RMSD_0'
  } else if(statS == 'MAE') {
    stat = ErrViewLib::mue
    ylab = paste0('MAE ',unit)
    if(normalize)
      ylab = 'MAE / MAE_0'
  } else {
    stat = ErrViewLib::q95hd
    ylab = paste0('Q95 ',unit)
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
    rep_probref  = rep_probref,
    plot         = plot,
    score        = score,
    confStat     = confStat,
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
#' Plot confidence curve for (uE,E) set. New version with improved validation.
#'
#' @param E (vector) prediction uncertainty, uE, or predicted value, V
#' @param uE (vector) error or z-score
#' @param probref (logical) plot probabilistic reference (probref) curve (default: FALSE)
#' @param rep_probref (integer) sampling repetitions for normal curve (default = 1000)
#' @param score (logical) estimate confidence score ? This requires a bootstrapped CI on the CC and might take some time.
#' @param nBoot (integer) nb of bootstrap repeats to estimate CI on confidence curve (default = 5000)
#' @param parallel (logical) use parallelized bootstrap ? (default = FALSE)
#' @param cl (object) cluster for parallel computing, which is used when parallel = TRUE. If parallel = TRUE and cl = NULL, then the cluster is defined as cl = makeCluster(detectCores()). (default = NULL)
#' @param plot (function) produce plot ?
#' @param col (integer) color index for the curve (default: 6)
#' @param add (logical) add confidence curve to previous plot (default: FALSE)
#' @param xlab (string) x axis label (default: 'k \% discarded)
#' @param xlim (vector) limits of the x axis (default: NULL)
#' @param ylim (vector) limits of the y axis (default: NULL)
#' @param title (string) a title to display above the plot (default: '')
#' @param showUk (logical) plot secondary axis with u_k values; supersedes `title`
#' @param showLegend (logical) display legend (default: TRUE)
#' @param legend (string) legend for the dataset (default: NULL)
#' @param unit (string) unit string to add to ylab
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
#'   plotCCVal(E, uE)
#' }
plotCCVal = function(
  E, uE,
  stat   = ErrViewLib::rmse,
  probref = FALSE,
  rep_probref = 1000,
  showUk = FALSE,
  plot   = TRUE,
  score  = FALSE,
  nBoot  = 5000,
  parallel = FALSE,
  cl     = NULL,
  col    = 5,
  type   = c('l','p'),
  add    = FALSE,
  xlab   = 'k% discarded',
  xlim   = NULL,
  ylab   = 'RMSE',
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

  if(score)
    probref = TRUE

  # Reorder data
  io = order(uE, decreasing = TRUE)
  uE = uE[io]
  E  = E[io]

  M = length(uE)
  pcVec = 0:99 # Vector of percentages

  # Statistics for data
  vstat = ErrViewLib::Sconf(E, pcVec, stat, FALSE)

  # Probabilistic reference   vnorm_mu = NULL
  if(probref) {
    # Simulated reference
    nrun = rep_probref
    vnorm = matrix(0, ncol = nrun, nrow = length(pcVec))
    for (i in 1:nrun) {
      X = uE * rnorm(M)
      vnorm[, i] = ErrViewLib::Sconf(X, pcVec, stat, FALSE)
    }
    vnorm_mu = apply(vnorm, 1, mean, na.rm = TRUE)
  }

  CS = CS_CI = NULL # Consistency score
  vstatCI = matrix(NA, nrow = length(pcVec), ncol = 2)
  if(score){
    ## Bootstrapped CI
    cl0 = NULL
    if(parallel) {
      library(parallel)
      cl0 = cl
      if(is.null(cl))
        cl <- makeCluster(detectCores())
    }
    for (i in 1:length(pcVec)) {
      k = pcVec[i]
      sel = 1:floor(k * M / 100)
      if (length(sel) == 0) {
        vstatCI[i,1:2] = NA
      } else {
        X0 = E[-sel]
        bs = nptest::np.boot(
          x = X0,
          statistic = stat, R = nBoot,
          level = 0.95, method = 'bca',
          parallel = parallel, cl = cl)
        vstatCI[i,1:2] = bs$bca
      }
    }
    if(parallel & is.null(cl0))
      stopCluster(cl)

    # Fraction of simulated CC within data CI
    success = sum(vnorm_mu >= vstatCI[,1] & vnorm_mu <= vstatCI[,2])
    trials  = length(pcVec)
    CS_CI   = DescTools::BinomCI(success, trials , method = "wilsoncc")
    CS      = CS_CI[1]
    ciPrint = paste0('[', round(CS_CI[2],2),', ',
                     round(CS_CI[3],2), ']')
  }

  if(plot) {
    # Expose gPars list
    for (n in names(gPars))
      assign(n, rlist::list.extract(gPars, n))

    if(add) {

      if(probref)
        lines(pcVec, vnorm_mu, lty = 3, lwd = 2*lwd, col=cols[col])

      if(score)
        polygon(
          c(pcVec,rev(pcVec)),
          c(vstatCI[,1],rev(vstatCI[,2])),
          col = cols_tr[col],
          border = NA)

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

      } else {
        title(main = title)

      }
      grid()

      if(probref)
        lines(pcVec, vnorm_mu, lty = 3, lwd = 2*lwd, col=cols[col])

      if(score)
        polygon(
          c(pcVec,rev(pcVec)),
          c(vstatCI[,1],rev(vstatCI[,2])),
          col = cols_tr[col],
          border = NA)

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

        if(probref) {
          legend = c('Prob. ref.',legend)
          lty    = c(3,lty)
          pch    = c(NA,pch)
          lcol   = c(col,lcol)
        }

        if(score) {
          legend = c(
            legend,
            paste0('CS = ',round(CS,2)),
            paste0('CI = ', ciPrint)
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

  }

  invisible(
    list(
      CS      = CS,
      CS_CI   = CS_CI,
      pcVec   = pcVec,
      curve   = vstat,
      curveCI = vstatCI,
      probRef = vnorm_mu
    )
  )

}
