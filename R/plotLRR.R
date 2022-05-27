#' Plot local range ratios
#'
#' @param E (vector) Errors
#' @param U (matrix) Prediction uncertainties
#' @param prob (vector) a set of coverage probabilities for the PUs
#' @param mycols (vector) a set of color indexes to gPars colors
#' @param nBin (integer) number of intervals for local coverage stats
#' @param ylim (vector) limits of the y axis
#' @param title (string) a title to display above the plot
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#' @param plot (logical) plot the results
#' @param ordX  (vector) set of abscissas to order sample
#' @param logX  (logical) log-transform abscissas
#' @param slide (logical) use sliding window
#' @param xlab  (string) abscissa label
#' @param xlim  (vector) range for abscissa
#' @param legend (string) legend for the dataset
#' @param legLoc (string) location of legend (see \link[grDevices]{xy.coord})
#' @param legNcol (integer) number of columns for legend
#' @param add (logical) add to previous graph ?
#'
#' @return Invisibly returns a list of LCP results. Mainly used
#'   for its plotting side effect.
#'
#' @export
#'
#' @examples
plotLRR = function(
  E, U,
  ordX      = NULL,
  logX      = FALSE,
  prob      = c(0.95),
  nBin      = NULL,
  plot      = TRUE,
  slide     = NULL,
  mycols    = 1:length(prob),
  xlab      = 'Calculated value',
  xlim      = NULL,
  ylim      = NULL,
  title     = '',
  legend    = paste0('P',round(100*prob)),
  legLoc    = 'bottom',
  legNcol   = 3,
  add       = FALSE,
  label     = 0,
  gPars     = ErrViewLib::setgPars()
) {


  N = length(E)
  if(NCOL(U) > 1) {
    if(length(prob) != NCOL(U))
      stop('>>> Provide target probabilities for all columns of U...')
  } else {
    U = as.matrix(U,ncol=1)
  }

  if(is.null(nBin))
    nBin  = max(min(floor(N/150),15),2)
  if(nBin <= 0)
    stop('>>> nBin should be > 0')
  if(is.null(slide))
    slide = nBin <= 4

  if(!is.null(ordX)) {
    if(length(ordX) != length(E))
      stop('>>> Inconsistent length for ordX')
    ord = order(ordX)
    xOrd = ordX[ord]
  } else {
    stop('>>> Please provide ordX')
  }
  eOrd = E[ord]
  uOrd = matrix(U[ord,],ncol = NCOL(U))

  # Design local areas
  intrv = ErrViewLib::genIntervals(N, nBin, slide)

  # Local Coverage stats
  pP = loP = upP = matrix(NA,nrow=length(prob),ncol=intrv$nbr)
  mint = c()

  # Width of p% empirical interval
  fRR = function(X, ind = 1:NROW(X), p = 0.95) {
    empi = diff(ErrViewLib::vhd(X[ind, 1], p = c((1 - p) / 2, (1 + p) / 2)))
    theo = 2 * sqrt(mean(X[ind, 2] ^ 2))
    return(theo / empi)
  }

  for (i in 1:intrv$nbr) {
    sel = intrv$lwindx[i]:intrv$upindx[i]
    err = eOrd[sel]
    M = length(sel)
    for (ip in seq_along(prob)) {
      p = prob[ip]
      unc = uOrd[sel,ip]
      bst = boot::boot(cbind(err,unc), fRR, R = 1000, p = p)
      bci = boot::boot.ci(bst, conf = 0.95, type = 'bca' )
      ci  = bci$bca[1, 4:5]
      pP[ip,i]  = bst$t0
      loP[ip,i] = ci[1]
      upP[ip,i] = ci[2]
    }
    mint[i] = 0.5*sum(range(xOrd[sel])) # Center of interval
  }

  if(plot) {
    # Plot ----
    if(length(gPars) == 0)
      gPars = ErrViewLib::setgPars()

    for (n in names(gPars))
      assign(n, rlist::list.extract(gPars, n))

    if(!add) {
      par(
        mfrow = c(1, 1),
        mar = mar, #c(mar[1:3],3),
        mgp = mgp,
        pty = 's',
        tcl = tcl,
        cex = cex,
        lwd = lwd
      )

      if(is.null(xlim))
        xlim = range(xOrd)

      if(is.null(ylim))
        ylim = range(c(loP, upP))

      matplot(
        mint,
        t(pP),
        xlab = xlab,
        ylab = 'Local Range Ratio',
        xlim = xlim,
        xaxs = 'i',
        ylim = ylim,
        yaxs = 'i',
        type = 'b',
        lty = 3,
        pch = 16,
        lwd = lwd,
        cex = ifelse(slide,0.5,1),
        col  = cols[mycols],
        main = title,
        log = ifelse(logX,'x','')
      )
      grid(equilogs = FALSE)
      abline(h   = 1,
             lty = 2,
             col = cols[1],
             lwd = lwd)

    } else {
      lines(
        mint,
        t(pP),
        type = 'b',
        lty = 3,
        pch = 16,
        lwd = lwd,
        cex = ifelse(slide,0.5,1),
        col  = cols[mycols]
      )
    }

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

    if(!add) {
      box()
      legend(
        legLoc, bty = 'n',
        legend = legend,
        col  = cols[mycols],
        lty  = 1,
        pch  = 16,
        ncol = legNcol,
        cex  = 0.8,
        adj  = 0.2
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
      mint   = mint,
      lwindx = intrv$lwindx,
      upindx = intrv$upindx,
      pc     = pP
      # pcl    = loP,
      # pcu    = upP,
      # meanP  = meanP,
      # uMeanP = uMeanP,
      # prob   = prob
    )
  )
}