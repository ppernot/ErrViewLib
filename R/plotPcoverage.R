#' Auxillary function for plotPcoverage
#'
#' @param cLearn
#' @param rLearn
#' @param cTest
#' @param rTest
#' @param corTrend
#' @param fo
#' @param prob
#' @param CImeth
#'
#' @return
#'
#' @examples
predIp = function(
  Learn, Test,
  corTrend = FALSE,
  fo       = NA,
  prob     = c(0.5, 0.75, 0.95),
  CImeth   = c('eq','pred','dist'),
  dist     = c('norm','t'),
  shape    = 2
){

  dist = match.arg(dist)
  CImeth = match.arg(CImeth)

  # Percentages for symmetric CI at prob values
  alpha = 1 - prob
  plow = alpha / 2
  psup = 1 - alpha / 2

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

  eqLwr = eqUpr = matrix(NA,nrow=length(prob),ncol=nrow(Test))

  if(CImeth == 'eq') {
    # Quantile estimates of CI limits
    for (k in seq_along(prob)) {
      eqLwr[k, ] = ErrViewLib::hd(E, q = plow[k])
      eqUpr[k, ] = ErrViewLib::hd(E, q = psup[k])
    }

  } else if(CImeth == 'pred') {
    # Linear regression of errors
    if (corTrend) {
      # Account for regression uncertainty
      for(k in seq_along(prob)) {
        cp = predict(reg,
                     newdata = Test,
                     interval = "prediction",
                     level = prob[k])
        eqUpr[k, ] = cp[,"upr"] - eP
        eqLwr[k, ] = cp[,"lwr"] - eP
      }
    } else {
      stop('CImeth = pred requires corTrend = TRUE')
    }
  } else {
    # Enlargement factor for unit variance distributions
    pfac = switch(
      dist,
      norm = normalp::qnormp(psup,p=shape) /
        sqrt(shape^(2/shape)*gamma(3/shape)/gamma(1/shape)),
      t    = qt(psup,df=shape) /
        sqrt(shape / (shape-2))
    ) # TBD: allow for non-symmetric distribs ?

    # Uniform CI from pdf hypothesis
    mu  = mean(E)
    sig = sd(E)
    for(k in seq_along(prob)) {
      eqLwr[k, ] = mu - pfac[k]*sig
      eqUpr[k, ] = mu + pfac[k]*sig
    }
  }

  return(
    list(
      eTest = eT,
      eqLwr = eqLwr,
      eqUpr = eqUpr
    )
  )
}
#' Plot local coverage probabilities to assess calibration and sharpness
#'
#' @param Data (data.frame) dataframe with predictor(s) and reference
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
#' @param nBin (integer) number of intervals for local coverage stats
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#' @param ordX  (vector) set of abscissas to order sample
#' @param logX  (logical) log-transform abscissas
#' @param binomCI (string) name of method to estimate Binomial Proportion CI
#' @param slide (logical) use sliding window
#' @param xlab  (string) abscissa label
#' @param xlim  (vector) range for abscissa
#' @param legLoc (string) location of legend (see \link[grDevices]{xy.coord})
#' @param legNcol (integer) number of columns for legend
#'
#' @return Invisibly returns a list of LCP results. Mainly used
#'   for its plotting side effect.
#'
#' @export
#'
#' @examples
plotPcoverage = function(
  # Reminder: changes to parameters should be reflected to alias plotLCP
  Data,
  corTrend  = FALSE,
  fo        = NA,
  ordX      = NULL,
  logX      = FALSE,
  CImeth    = c('eq','pred','dist'),
  prob      = c(0.5,0.75,0.95),
  dist      = c('norm','t'),
  shape     = 2,
  valid     = c("kfold","loo"),
  nFold     = 10,
  nRepeat   = 10,
  nBin      = NULL,
  binomCI   = c("wilson", "wilsoncc", "clopper-pearson",
               "agresti-coull", "jeffreys"),
  plot      = TRUE,
  slide     = NULL,
  mycols    = 1:length(prob),
  xlab      = 'Calculated value',
  xlim      = NULL,
  ylim      = c(0,1),
  title     = '',
  legLoc    = 'bottom',
  legNcol   = 3,
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {

  dist    = match.arg(dist)
  valid   = match.arg(valid)
  CImeth  = match.arg(CImeth)
  binomCI = match.arg(binomCI)
  # ordX    = match.arg(ordX)

  if(nRepeat <= 0)
    stop('>>> nRepeat should be > 0')
  if(nFold <= 0)
    stop('>>> nFold should be > 0')


  C = Data$D
  R = Data$R
  E = R - C
  Data = cbind(Data,E)
  N = length(C)

  if(is.null(nBin))
    nBin  = max(min(floor(N/150),15),2)
  if(nBin <= 0)
    stop('>>> nBin should be > 0')
  if(is.null(slide))
    slide = nBin <= 4

  ord = order(C)
  xOrd = C[ord]
  # if( ordX == "uP" & !is.null(Data$uP)) {
  #   ord = order(Data$uP)
  #   xOrd = Data$uP[ord]
  # }
  if(!is.null(ordX)) {
    if(length(ordX) != length(C))
      stop('>>> Inconsistent length for ordX')
    ord = order(ordX)
    xOrd = ordX[ord]
  }
  DataOrd = Data[ord,]

  # Design local areas
  intrv  = ErrViewLib::genIntervals(N, nBin, slide)
  nbr    = intrv$nbr
  lwindx = intrv$lwindx
  upindx = intrv$upindx

  if(!is.null(Data$uP)) {
    # z-scores ----
    uP = Data$uP
    zOrd = (R[ord]-C[ord])/uP[ord]

    # Direct validation of z-scores: no cross-validation

    alpha = 1 - prob
    plow  = alpha / 2
    psup  = 1 - alpha / 2
    qlow  = switch(
      dist,
      norm = normalp::qnormp(plow,p=shape) /
        sqrt(shape^(2/shape)*gamma(3/shape)/gamma(1/shape)),
      t    = qt(plow,df=shape) /
        sqrt(shape / (shape-2))

    )
    qsup  = switch(
      dist,
      norm = normalp::qnormp(psup,p=shape) /
        sqrt(shape^(2/shape)*gamma(3/shape)/gamma(1/shape)),
      t    = qt(psup,df=shape) /
        sqrt(shape / (shape-2))

    )
    # Rq: upper and lower limits are maintained for
    # future extension to skewed distribution...

    # Generate coverage tests matrix
    tm = matrix(NA,ncol=length(prob),nrow=length(zOrd))
    for (ip in seq_along(prob))
      tm[,ip] = as.numeric(zOrd >= qlow[ip] & zOrd <= qsup[ip])

    nRepeat = 1

  } else {
    # Errors ----

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

        Learn = DataOrd[iLearn,]
        Test  = DataOrd[iTest,]

        # Prediction of CI over eTest
        predCI = predIp(
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
              predCI$eTest[j] >= predCI$eqLwr[ip, j] &
              predCI$eTest[j] <= predCI$eqUpr[ip, j]
            )
            # Accumulate for kfold repeats
            tm[iTest[j], ip] = tm[iTest[j], ip] + t
          }
        }
      } # End nFold loop: all points tested
    } # End nRepeat loop
  } # End z-score test

  # Coverage stats for test matrix
  pP = loP = upP = matrix(NA,nrow=length(prob),ncol=nbr)
  mint = c()
  for (i in 1:nbr) {
    sel = lwindx[i]:upindx[i]
    M = length(sel)
    for (ip in seq_along(prob)) {
      S = sum(tm[sel,ip])/nRepeat
      pP[ip,i] = S / M
      ci = DescTools::BinomCI(S, M, method = binomCI)
      loP[ip, i] = ci[,2]
      upP[ip, i] = ci[,3]
    }
    mint[i] = 0.5*sum(range(xOrd[sel])) # Center of interval

  }
  S = colSums(tm)/nRepeat
  ci = DescTools::BinomCI(S, N, method = binomCI)
  meanP = ci[,1]
  uMeanP = sqrt(meanP*(1-meanP)/N)
  loMeanP = ci[,2]
  upMeanP = ci[,3]

  if(plot) {
    # Plot ----
    if(length(gPars) == 0)
      gPars = ErrViewLib::setgPars()

    for (n in names(gPars))
      assign(n, rlist::list.extract(gPars, n))

    par(
      mfrow = c(1, 1),
      mar = c(mar[1:3],3),
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd
    )

    if(is.null(xlim))
      xlim = range(xOrd)

    if (any(is.na(ylim)))
      ylim = range(c(loP, upP))

    matplot(
      mint,
      t(pP),
      xlab = xlab,
      ylab = 'Local Coverage Probability',
      xlim = xlim,
      xaxs = 'i',
      ylim = ylim,
      type = 'b',
      lty = 3,
      pch = 16,
      lwd = lwd,
      cex = ifelse(slide,0.5,1),
      col  = cols[mycols],
      main = title,
      log = ifelse(logX,'x','')
    )
    grid()

    if(slide) {
      ipl = seq(1,length(mint),length.out=nBin)
      for(i in seq_along(prob)) {
        polygon(
          c(mint,rev(mint)),
          c(loP[i,], rev(upP[i,])),
          col = cols_tr[mycols[i]],
          border = NA)
        segments(
          mint[ipl], loP[i,ipl],
          mint[ipl], upP[i,ipl],
          col  = cols[mycols[i]],
          lwd  = 1.5 * lwd,
          lend = 1)
      }

    } else {
      for(i in seq_along(prob))
        segments(
          mint, loP[i,],
          mint, upP[i,],
          col  = cols[mycols[i]],
          lwd  = 1.5 * lwd,
          lend = 1)

    }
    mtext(text = paste0(prob,' -'),
          side = 2,
          at = prob,
          col = cols[mycols],
          cex = 0.75*cex,
          las = 1,
          font = 2)
    xpos = pretty(xOrd)
    abline(h   = prob,
           lty = 2,
           col = cols[mycols],
           lwd = lwd)

    box()

    # Mean coverage proba
    ypos = par("usr")[4]
    pm = c()
    for(i in seq_along(meanP)) {
      if(uMeanP[i] != 0)
        pm[i] = ErrViewLib::prettyUnc(meanP[i],uMeanP[i],1)
      else
        pm[i] = signif(meanP[i],2)
    }
    mtext(text = c(' Mean',paste0('- ',pm)),
          side = 4,
          at = c(ypos,meanP),
          col = c(1,cols[mycols]),
          cex = 0.75*cex,
          las = 1,
          font = 2)
    for(i in seq_along(meanP))
      segments(
        xlim[2],loMeanP[i],
        xlim[2],upMeanP[i],
        col  = cols[mycols][i],
        lwd  = 6 * lwd,
        lend = 1
      )

    legend(
      legLoc, bty = 'n',
      legend = paste0('P',round(100*prob)),
      col  = cols[mycols],
      lty  = 1,
      pch  = 16,
      ncol = legNcol,
      cex  = 0.8,
      adj  = 0.2
    )

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
      mint   = mint,
      lwindx = lwindx,
      upindx = upindx,
      pc     = pP,
      pcl    = loP,
      pcu    = upP,
      meanP  = meanP,
      uMeanP = uMeanP,
      prob   = prob
    )
  )
}