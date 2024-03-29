#' Auxillary function for plotCalCurve
#' Predicted quantiles of errors
#'
#' @param cLearn -
#' @param rLearn -
#' @param cTest -
#' @param rTest -
#' @param corTrend -
#' @param fo -
#' @param prob -
#' @param CImeth -
#'
#' @return a list
#'
predQ = function(
  Learn, Test,
  corTrend = FALSE,
  fo       = NA,
  prob     = seq(0.005, 0.995, by=0.005),
  CImeth   = c('eq','pred','dist'),
  dist     = c('norm','t'),
  shape    = 2
){

  dist = match.arg(dist)
  CImeth = match.arg(CImeth)

  # Errors
  E  = Learn[,'R'] - Learn[,'D']
  eT = Test[,'R']  - Test[,'D']

  # Apply trend corrections
  if(corTrend) {
    environment(fo) <- environment()
    reg = lm(fo, data = Learn)
    E  = residuals(reg)
    eP  = predict(reg, newdata = Test)
    eT  = eT - eP
  }

  eqLwr = matrix(NA,nrow=length(prob),ncol=nrow(Test))

  if(CImeth == 'eq') {
    # Quantile estimates of CI limits
    for (k in seq_along(prob))
      eqLwr[k, ] = ErrViewLib::hd(E, q = prob[k])

  } else if(CImeth == 'pred') {
    # Linear regression of errors
    if (corTrend) {
      # Account for regression uncertainty
      for(k in seq_along(prob)) {
        # df = ciTools::add_quantile(
        #   Test,
        #   reg,
        #   p=prob[k],
        #   name='Q')
        #   eqLwr[k, ] = df[,'Q'] - eP
        #   # Avoid ciTools because of build problems in docker
        out <- predict(
          reg, newdata = Test,
          interval = "prediction",
          se.fit = TRUE)
        fitted <- out$fit[,1]
        residual_df <- out$df
        se_fitted <- out$se.fit
        resid_var <- out$residual.scale^2
        se_pred <- sqrt(resid_var + se_fitted^2)
        t_quantile <- qt(p = prob[k], df = residual_df)
        out_quantiles <- fitted + se_pred * t_quantile
        eqLwr[k, ] = out_quantiles - eP
      }
    } else {
      stop('"CImeth = pred" requires "corTrend = TRUE"')
    }
  } else {
    pfac = switch(
      dist,
      norm = normalp::qnormp(prob,p=shape) /
        sqrt(shape^(2/shape)*gamma(3/shape)/gamma(1/shape)),
      t    = qt(prob,df=shape) /
        sqrt(shape / (shape-2))
    )
    # TBD: allow for non-symmetric distribs ?

    # Quantiles from pdf hypothesis
    mu  = mean(E)
    sig = sd(E)
    for(k in seq_along(prob))
      eqLwr[k, ] = mu + pfac[k]*sig

  }

  return(
    list(
      eTest = eT,
      eqLwr = eqLwr
    )
  )
}
#' Older version of calibration curves. Kept for reproducibility of papers.
#'
#' @param Data (data.frame) dataframe with predictor(s) and reference
#' @param uP (vector) a set of prediction uncertainties
#' @param corTrend (logical) flag to correct trend
#' @param fo (formula) trend correction formula
#' @param CImeth (string) method to estimate CI limits
#' @param prob (vector) a set of coverage probabilities to test
#' @param dist (string) a distribution
#' @param shape (numeric) shape parameter (> 0, maybe non-integer).
#' @param mycols (vector) a set of color indexes to gPars colors
#' @param valid (string) a cross-validation method: "kfold" (default) or "loo"
#' @param nFold (integer) number of folds (default: 10)
#' @param nRepeat (integer) number of repeats for k-fold cross-validation
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#'
#' @return Invisibly returns a list of calibration statistics. Mainly used
#'   for its plotting side effect.
#' @export
#'
plotCalCurve = function(
  Data,
  corTrend  = FALSE,
  fo        = NA,
  CImeth    = c('eq','pred','dist'),
  prob      = seq(0.005, 0.995, by=0.005),
  dist      = c('norm','t'),
  shape     = 2,
  valid     = c("kfold","loo"),
  nFold     = 10,
  nRepeat   = 10,
  score     = TRUE,
  binomCI   = c("wilson", "wilsoncc", "clopper-pearson",
               "agresti-coull", "jeffreys"),
  title     = '',
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  dist    = match.arg(dist)
  valid   = match.arg(valid)
  CImeth  = match.arg(CImeth)
  binomCI = match.arg(binomCI)

  if(nRepeat <= 0)
    stop('>>> nRepeat should be > 0')
  if(nFold <= 0)
    stop('>>> nFold should be > 0')


  C = Data$D
  R = Data$R
  E = R - C
  Data = cbind(Data,E)
  N = length(C)


  if(valid =="loo") {
    cfold = 1:N
    nFold = N
    nRepeat = 1

  } else {
    # Generate kfold sample of approximately same size
    # enables to distribute randomly points above round
    # division.
    pfold = seq(0, 1, length.out = nFold + 1)[-1]
    bfold = vhd(1:N, p = pfold)
    cfold = c()
    for (i in 1:N)
      cfold[i] = which(bfold >= i)[1]
    if (sum(table(cfold)) != N)
      stop('>>> Pb. in cfold calculation ----')

  }

  # Generate samples
  tm  = matrix(0, nrow = N, ncol = length(prob))
  for(irep in 1:nRepeat) {
    iran = sample.int(N,N) # Randomize points
    for (k in 1:nFold) {
      sel    = which(cfold == k)
      iTest  = iran[sel]
      iLearn = which(! iran %in% iTest)

      Learn = Data[iLearn,]
      Test  = Data[iTest,]

      # Prediction of CI over eTest
      predCI = predQ(
        Learn, Test,
        corTrend = corTrend,
        fo = fo,
        prob = prob,
        CImeth = CImeth,
        dist = dist,
        shape = shape
      )

      # Test
      for (ip in seq_along(prob)) {
        for (j in 1:nrow(Test)) {
          t = as.numeric(
            predCI$eTest[j] <= predCI$eqLwr[ip, j]
          )
          # Accumulate for kfold repeats
          tm[iTest[j], ip] = tm[iTest[j], ip] + t
        }
      }
    } # End nFold loop: all points tested
  } # End nRepeat loop

  # Coverage stats for test matrix
  pP = c()
  for (ip in seq_along(prob))
    pP[ip] = sum(tm[,ip]) / nRepeat / N

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = 's',
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  pt = c(0,prob,1)
  pe = c(0,pP,1)
  plot(
    pt, pe,
    type = 'l',
    xlim = c(0,1),
    xlab = 'Expected cumulative distribution',
    ylim = c(0,1),
    ylab = 'Observed cumulative distribution',
    xaxs = 'i',
    yaxs = 'i',
    col  = cols[3],
    main = title,
    lwd  = lwd
  )
  grid()
  abline(
    a=0,b=1,
    lty=2,
    col = cols[6],
    lwd = lwd
  )
  polygon(
    c(pt,0),c(pe,0),
    col = cols_tr[5],
    border = NA
  )
  lwr = upr = c()
  for(i in seq_along(pe)) {
    ci = DescTools::BinomCI(pe[i]*N, N, method = 'wilsoncc')
    lwr[i] = ci[,2]
    upr[i] = ci[,3]
  }
  matlines(pt,cbind(lwr,upr),col = cols[3], lty = 2, lwd = lwd)
  box()
  if(score) {
    trapz = function(x,y) {
      # Trapezoidal integration
      idx = 2:length(x)
      return (as.double( (x[idx] - x[idx-1]) %*%
                           (y[idx] + y[idx-1])) / 2)
    }
    misCal   = trapz(pt,abs(pe-pt))
    misCalUp = trapz(pt,abs(upr-lwr)) / 2
    calErr   = sqrt(sum((pt-pe)^2))
    text(
      0.75, 0.1,
      paste0(
        'MisCal   = ',signif(misCal,2),'\n',
        'MisCalUp = ',signif(misCalUp,2),'\n',
        'CalErr   = ',signif(calErr,2)
      ),
      cex=1
    )
  }

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)

  invisible(
    list(
      misCal = misCal,
      misCalUp = misCalUp,
      calErr = calErr
    )
  )
}
#' Plot calibration curve
#'
#' @param E (vector) a vector of errors
#' @param uE (vector) a vector of uncertainties
#' @param dist (string) a distribution (default: `Normal`)
#' @param shape (numeric) shape parameter for the T and Normp distributions
#' @param score (logical) evaluate the calibration scores (default: `TRUE`)
#' @param plot (logical) plot the results (default: `TRUE`)
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#'
#' @return Invisibly returns a list of calibration statistics. Mainly used
#'   for its plotting side effect.
#' @export
#'
plotCalibration = function(
  E, uE,
  prob      = seq(0.005, 0.995, by=0.005),
  dist      = c('Normal','Student','Uniform','Laplace','Normp','T4','Normp4'),
  shape     = 4,
  score     = TRUE,
  plot      = TRUE,
  title     = '',
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  dist    = match.arg(dist)

  N = length(E)

  # Coverage stats data
  Normal = function(p, ...)
    qnorm(p)
  Student = function(p, df = 4)
    qt(p, df = df) / sqrt(df/(df-2))
  Uniform = function(p, ...)
    qunif(p, -sqrt(3), sqrt(3))
  Laplace = function(pr, ...)
    normalp::qnormp(pr, p = 1) / sqrt(gamma(3)/gamma(1))
  Normp = function(pr, df = 4)
    normalp::qnormp(pr, p = df) / sqrt(df^(2/df)*gamma(3/df)/gamma(1/df))
  T4 = function(p, ...)
    Student(p, df = 4)
  Normp4 = function(p, ...)
    Normp(p, df = 4)

  qfun = get(dist)

  pP = c()
  for (ip in seq_along(prob))
    pP[ip] = mean(E <= uE * qfun(prob[ip], df = shape))

  pt = c(0,prob,1)
  pe = c(0,pP,1)

  misCal = misCalUp = calErr = NA
  if(score) {
    trapz = function(x,y) {
      # Trapezoidal integration
      idx = 2:length(x)
      return (as.double( (x[idx] - x[idx-1]) %*%
                           (y[idx] + y[idx-1])) / 2)
    }
    lwr = upr = c()
    for(i in seq_along(pe)) {
      ci = DescTools::BinomCI(pe[i]*N, N, method = 'wilsoncc')
      lwr[i] = ci[,2]
      upr[i] = ci[,3]
    }
    misCal   = trapz(pt,abs(pe-pt))
    misCalUp = trapz(pt,abs(upr-lwr)) / 2
    calErr   = sqrt(sum((pt-pe)^2))
  }

  if(plot) {

    # Expose gPars list
    for (n in names(gPars))
      assign(n, rlist::list.extract(gPars, n))

    par(
      mfrow = c(1, 1),
      mar = mar,
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd
    )

    plot(
      pt, pe,
      type = 'l',
      xlim = c(0,1),
      xlab = 'Expected cumulative distribution',
      ylim = c(0,1),
      ylab = 'Observed cumulative distribution',
      xaxs = 'i',
      yaxs = 'i',
      col  = cols[3],
      main = title,
      lwd  = lwd
    )
    grid()
    abline(
      a=0,b=1,
      lty=2,
      col = cols[6],
      lwd = lwd
    )
    polygon(
      c(pt,0),c(pe,0),
      col = cols_tr[5],
      border = NA
    )
    lwr = upr = c()
    for(i in seq_along(pe)) {
      ci = DescTools::BinomCI(pe[i]*N, N, method = 'wilsoncc')
      lwr[i] = ci[,2]
      upr[i] = ci[,3]
    }
    matlines(pt,cbind(lwr,upr),col = cols[3], lty = 2, lwd = lwd)
    box()
    if(score) {
      text(
        0.75, 0.1,
        paste0(
          'MisCal   = ',signif(misCal,2),'\n',
          'MisCalUp = ',signif(misCalUp,2),'\n',
          'CalErr   = ',signif(calErr,2)
        ),
        cex=1
      )
    }

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
      misCal = misCal,
      misCalUp = misCalUp,
      calErr = calErr
    )
  )
}
#' Plot calibration curve from quantiles table
#'
#' @param E (vector) a vector of errors
#' @param Q (list) a list of two elements: `Q$prob` is a vector of probabilities #'    and `Q$quant` is a matrix of quantile values for each probability (column) #{    and uncertainty (row)
#' @param cond (vector) a vector of values on which conditional calibration curves are designed
#' @param intrv (object) intervals generated by `genIntervals` (default: `NULL`)
#' @param nBin (integer) number of intervals for local coverage stats
#' @param equiPop (logical) generate intervals  with equal bin counts
#'   (default: `equiPop = TRUE`)
#' @param popMin (integer) minimal bin count in an interval
#' @param logBin (logical) if `equiPop = FALSE`, one can choose between
#'   equal range intervals, or equal log-range intervals
#'   (default `logBin = TRUE`)
#' @param score (logical) evaluate the calibration scores (default: `TRUE`)
#' @param plot (logical) plot the results (default: `TRUE`)
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#'
#' @return Invisibly returns a list of calibration statistics. Mainly used
#'   for its plotting side effect.
#' @export
#'
plotCalibrationQ = function(
  E, Q,
  cond      = NULL,
  nBin      = NULL,
  intrv     = NULL,
  equiPop   = TRUE,
  popMin    = 30,
  logBin    = TRUE,
  score     = TRUE,
  plot      = TRUE,
  title     = '',
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  N = length(E)

  quant = Q$quant
  prob  = Q$prob

  # Conditional cal-curves
  plotCond = !is.null(cond) & is.vector(cond)
  if(plotCond) {
    iord = order(cond)
    xOrd = cond[iord]
    # Design local areas
    if(is.null(intrv)) {
      if(is.null(nBin))
        nBin  = max(min(floor(N/150),15),2)
      if(nBin <= 0)
        stop('>>> nBin should be > 0')
      Xin = N
      if(!equiPop)
        Xin = xOrd
      slide = FALSE
      intrv = ErrViewLib::genIntervals(Xin, nBin, slide, equiPop,
                                       popMin, logBin)
    }
    nBin0 = nBin
    nBin = intrv$nbr
    peCond = matrix(NA, nrow = length(prob) + 2, ncol = nBin)
    for (i in 1:nBin) {
      sel = intrv$lwindx[i]:intrv$upindx[i]
      pP = c()
      for (ip in seq_along(prob))
        pP[ip] = mean(E[iord][sel] <= quant[iord,ip][sel])
      peCond[,i] = c(0,pP,1)
    }
  }

  # Average cal-curve
  pP = c()
  for (ip in seq_along(prob))
    pP[ip] = mean(E <= quant[,ip])

  # Add extreme values
  pt = c(0,prob,1)
  pe = c(0,pP,1)

  # Calibration statistics
  misCal = misCalUp = calErr = NA
  if(score) {
    trapz = function(x,y) {
      # Trapezoidal integration
      idx = 2:length(x)
      return (as.double( (x[idx] - x[idx-1]) %*%
                           (y[idx] + y[idx-1])) / 2)
    }
    lwr = upr = c()
    for(i in seq_along(pe)) {
      ci = DescTools::BinomCI(pe[i]*N, N, method = 'wilsoncc')
      lwr[i] = ci[,2]
      upr[i] = ci[,3]
    }
    misCal   = trapz(pt,abs(pe-pt))
    misCalUp = trapz(pt,abs(upr-lwr)) / 2
    calErr   = sqrt(sum((pt-pe)^2))
  }

  if(plot) {

    # Expose gPars list
    for (n in names(gPars))
      assign(n, rlist::list.extract(gPars, n))

    par(
      mfrow = c(1, 1),
      mar = mar,
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd
    )

    plot(
      pt, pe,
      type = 'l',
      xlim = c(0,1),
      xlab = 'Expected cumulative distribution',
      ylim = c(0,1),
      ylab = 'Observed cumulative distribution',
      xaxs = 'i',
      yaxs = 'i',
      col  = cols[5],
      main = title,
      lwd  = 2*lwd
    )
    grid()
    nCond = 1
    if(plotCond) {
      matplot(
        pt, peCond,
        type = 'l',
        col  = cols_tr2[5],
        lty  = 1,
        lwd  = 1.5*lwd,
        add = TRUE
      )
      nCond = nBin
    }
    abline(
      a=0,b=1,
      lty = 2,
      col = cols[1],
      lwd = 1.5*lwd
    )
    polygon(
      c(pt,0),c(pe,0),
      col = cols_tr[5],
      border = NA
    )
    # Confidence area
    lwr = upr = c()
    if(plotCond) {
      for(i in seq_along(pt)) {
        ci = DescTools::BinomCI(pt[i]*N/nCond, N/nCond, method = 'wilsoncc')
        lwr[i] = ci[,2]
        upr[i] = ci[,3]
      }
    } else {
      for(i in seq_along(pe)) {
        ci = DescTools::BinomCI(pe[i]*N, N, method = 'wilsoncc')
        lwr[i] = ci[,2]
        upr[i] = ci[,3]
      }
    }
    matlines(pt,cbind(lwr,upr),col = cols[3], lty = 2, lwd = 1.5*lwd)
    box()
    if(score) {
      text(
        0.75, 0.1,
        paste0(
          'MisCal   = ',signif(misCal,2),'\n',
          'MisCalUp = ',signif(misCalUp,2),'\n',
          'CalErr   = ',signif(calErr,2)
        ),
        cex=1
      )
    }

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
      misCal = misCal,
      misCalUp = misCalUp,
      calErr = calErr
    )
  )
}
