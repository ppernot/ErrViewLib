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
predQ = function(
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
    pfac = switch(
      dist,
      norm = normalp::qnormp(psup,p=shape),
      t    = qt(psup,df=shape)
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
#' @param nBin (integer) number of intervals for local coverage stats
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param legloc (string) location of the legend (default: "bottom")
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#'
#' @return Invisibly returns a list of LCP results. Mainly used
#'   for its plotting side effect.
#' @export
#'
#' @examples
plotPcoverage = function(
  Data,
  corTrend  = FALSE,
  fo        = NA,
  CImeth    = c('eq','pred','dist'),
  prob      = c(0.5,0.75,0.95),
  dist      = c('norm','t'),
  shape     = 2,
  valid     = c("kfold","loo"),
  nFold     = 10,
  nRepeat   = 10,
  nBin      = 10,
  plot      = TRUE,
  mycols    = 1:length(prob),
  xlab      = 'Calculated value',
  ylim      = NA,
  title     = '',
  legloc    = 'bottom',
  label     = 0,
  gPars     = NULL
) {

  dist  = match.arg(dist)
  valid = match.arg(valid)
  CImeth = match.arg(CImeth)

  if(nRepeat <= 0)
    stop('>>> nRepeat should be > 0')
  if(nFold <= 0)
    stop('>>> nFold should be > 0')
  if(nBin <= 0)
    stop('>>> nBin should be > 0')

  C = Data$D
  R = Data$R
  E = R - C
  Data = cbind(Data,E)
  N = length(C)

  if(!is.null(Data$uP)) {
    uP = Data$uP

    # z-scores ----
    # Direct validation of z-scores: no cross-validation

    alpha = 1 - prob
    plow  = alpha / 2
    psup  = 1 - alpha / 2
    qlow  = switch(
      dist,
      norm = normalp::qnormp(plow, p = shape),
      t    = qt(plow, df = shape)
    )
    qsup  = switch(
      dist,
      norm = normalp::qnormp(psup, p = shape),
      t    = qt(psup, df = shape)
    )

    # Attribute bin numbers to data
    ord  = order(C)
    cOrd = C[ord]
    zOrd = (R[ord]-C[ord])/uP[ord]
    p    = seq(0, 1, length.out = nBin + 1)[-1]
    br   = vhd(cOrd, p = p)
    cl   = c()
    for (i in seq_along(cOrd))
      cl[i] = which(br >= cOrd[i])[1]

    # Coverage stats
    pP = loP = upP = matrix(NA,nrow=length(prob),ncol=length(br))
    mint = c()
    tG = matrix(NA,nrow=length(prob),ncol=length(cl))
    i0 = 1
    for (i in seq_along(br)) {
      sel = which(cl==i)
      len = length(sel)
      for (ip in seq_along(prob)) {
        t = zOrd[sel] >= qlow[ip] & zOrd[sel] <= qsup[ip]
        tG[ip,i0:(i0+len-1)] = t
        pp         = mean(t)
        pP[ip,i]   = pp
        limits     = binom::binom.confint(
          x = sum(t),
          n = length(t),
          methods = "wilson")
        loP[ip, i] = limits$lower
        upP[ip, i] = limits$upper
      }
      i0 = i0 + len
      mint[i] = mean(cOrd[sel]) # Center of interval
    }
    meanP = rowMeans(tG) # avoid inequal samples bias
    cvP = apply(pP,1,sd)/meanP*100

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
    pin  = matrix(0, nrow = N, ncol = length(prob))
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
            t = predCI$eTest[j] >= predCI$eqLwr[ip, j] &
              predCI$eTest[j] <= predCI$eqUpr[ip, j]
            # Accumulate for kfold repeats
            pin[iTest[j], ip] = pin[iTest[j], ip] + t
          }
        }

      } # End nFold loop: all points tested

    } # End nRepeat loop

    # Partition of predictive variable for local coverage stats
    ord = order(C)
    cOrd = C[ord]
    p  = seq(0, 1, length.out = nBin + 1)[-1]
    br = vhd(cOrd, p = p)
    cl = c()
    for (i in seq_along(cOrd))
      cl[i] = which(br >= cOrd[i])[1]

    # Coverage stats
    pP = loP = upP = matrix(NA,nrow=length(prob),ncol=length(br))
    mint = c()
    for (i in seq_along(br)) {
      sel = which(cl==i)
      for (ip in seq_along(prob)) {
        X          = pin[ord[sel],ip]
        pp         = ErrViewLib::mse(X) / nRepeat
        pP[ip,i]   = pp
        limits     = binom::binom.confint(
          x = sum(X),
          n = length(X)*nRepeat,
          methods = "wilson")
        loP[ip, i] = limits$lower
        upP[ip, i] = limits$upper
      }
      mint[i] = mean(cOrd[sel]) # Center of interval
    }
    meanP = rowMeans(pP) # Mean coverage
    cvP = apply(pP,1,sd)/meanP*100

  }

  if(plot) {

    # Plot ----
    if(length(gPars) == 0)
      gPars = ErrViewLib::setgPars()

    for (n in names(gPars))
      assign(n, rlist::list.extract(gPars, n))

    par(
      mfrow = c(1, 1),
      mar = c(mar[1:3],4.3),
      mgp = mgp,
      pty = 's',
      tcl = tcl,
      cex = cex,
      lwd = lwd
    )

    if (any(is.na(ylim)))
      ylim = range(c(loP, upP))

    matplot(
      mint,
      t(pP),
      xlab = xlab,
      ylab = 'Local Coverage Probability',
      ylim = ylim,
      type = 'b',
      lty = 3,
      pch  = 19,
      lwd = lwd,
      col  = cols[mycols],
      main = title
    )
    grid()
    for(i in seq_along(prob))
      segments(mint, loP[i,], mint, upP[i,], col = cols[mycols[i]])
    mtext(text = paste0(prob,' -'),
          side = 2,
          at = prob,
          col = cols[mycols],
          cex = 0.75*cex,
          las = 1,
          font = 2)
    abline(h   = prob,
           lty = 2,
           col = cols[mycols],
           lwd = lwd)
    # Mean coverage proba
    mtext(text = c(' Mean',paste0('- ',signif(meanP,2))),
          side = 4,
          at = c(1.02,meanP),
          col = c(1,cols[mycols]),
          cex = 0.75*cex,
          las = 1,
          font = 2)
    mtext(text = c('(CV)',paste0('(',signif(cvP,2),' %)')),
          side = 4,
          at = c(1.02,meanP),
          col = c(1,cols[mycols]),
          cex = 0.75*cex,
          las = 1,
          line = 2,
          font = 2)
    legend(
      legloc, bty = 'n',
      legend = paste0('P',round(100*prob)),
      col = cols[mycols],
      pch = 19,
      lty  = 3,
      ncol = 1,
      cex  = 0.8
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

  invisible(
    list(
      pc    = pP,
      pcl   = loP,
      pcu   = upP,
      meanP = meanP,
      cvP   = cvP,
      prob  = prob
    )
  )
}