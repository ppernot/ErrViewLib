#' Plot of the CDF of the differences of a pair of absolute errors samples
#'
#'Plots also auxiliary statistics (SIP, MG, ML) with bootstrapped CIs.
#'
#' @param err
#' @param meth1
#' @param meth2
#' @param eps
#' @param xmax
#' @param xlab
#' @param units
#' @param main
#' @param nboot
#' @param label
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotDeltaCDF <- function(
  err,
  meth1,
  meth2,
  eps   = NULL,
  xmin  = NULL,
  xmax  = NULL,
  xlab  = NULL,
  units = 'a.u.',
  main  = '',
  nboot = 1000,
  label = 0,
  showSIP  = TRUE,
  showMLG  = TRUE,
  showDmue = TRUE,
  showCI   = TRUE,
  gPars
) {
  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2,
    xaxs = 'i',
    yaxs = 'i'
  )

  if(class(err)=='matrix' | class(err)=='data.frame') {
    X = abs(err[,meth1]) - abs(err[,meth2])
    Y = cbind(err[,meth1],err[,meth2])
  } else {
    X = abs(err[[meth1]]) - abs(err[[meth2]])
    Y = cbind(err[[meth1]],err[[meth2]])
  }

  nboot = max(nboot, length(X)+1)

  # Stats of SIP and DeltaMUE indicators
  if (showSIP | showMLG) {
    bsSip     = boot::boot(Y, statistic = fsi, R = nboot)
    if (showCI) {
      bsSip.ci  = boot::boot.ci(bsSip,
                                index = 1,
                                conf = 0.95,
                                type = "bca")
      bsMG.ci   = boot::boot.ci(bsSip,
                                index = 2,
                                conf = 0.95,
                                type = "bca")
      bsML.ci   = boot::boot.ci(bsSip,
                                index = 3,
                                conf = 0.95,
                                type = "bca")
    }
  }
  if (showDmue) {
    bsDmue    = boot::boot(Y, statistic = dmue, R = nboot)
    if (showCI) {
      bsDmue.ci = boot::boot.ci(bsDmue,
                                conf = 0.95,
                                type = "bca")
    }
  }

  if(is.null(xmin))
    xmin = min(X)
  if(is.null(xmax))
    xmax = max(X)
  xlim = c(xmin,xmax)

  if(is.null(xlab))
    xlab = substitute(
      abs(Error[meth1])-abs(Error[meth2])~~group("[",units,"]"),
      list(meth1 = meth1, meth2 = meth2, units = units))

  x = sort(X,na.last = NA)
  y = (1:length(x))/length(x)
  plot(
    x, y,
    type = 'l',
    col  = cols[5],
    xlab = xlab,
    xlim = xlim,
    ylab = 'Probability'
  )
  title(main=main,cex.main=0.9,adj=0,line=1)
  grid()
  if(showCI) {
    # Approx. CI on ECDF
    sigp = sqrt(y * (1 - y) / length(y))
    polygon(c(x, rev(x)),
            c(y - 1.96 * sigp, rev(y + 1.96 * sigp)),
            col = cols_tr2[5],
            border = NA)
  }
  if (!is.null(eps))
    polygon(
      x = c(-eps, -eps, eps, eps),
      y = c(0, 1, 1, 0),
      col = cols_tr2[3],
      border = NA
    )
  abline(v = 0,col=1,lty=1)

  # MG / ML
  if(showMLG) {
    mval = bsSip$t0
    if(showCI) {
      q025 = q975 = c()
      q025[2] = bsMG.ci$bca[,4]
      q975[2] = bsMG.ci$bca[,5]
      q025[3] = bsML.ci$bca[,4]
      q975[3] = bsML.ci$bca[,5]
      for(i in 2:3)
        polygon(
          c(q025[i],q975[i],q975[i],q025[i]),
          c(0,0,1,1),
          col = cols_tr2[4],
          border = NA)
    }
    abline(v = mval[2:3], col = cols[6], lty=2)
    mtext(text = c('MG','ML'),
          side = 3,
          at   = mval[2:3],
          col  = cols[6],
          cex  = 0.75*cex)

  }
  if(showSIP) {
    mval = bsSip$t0
    if(showCI) {
      q025 = q975 = c()
      q025[1] = bsSip.ci$bca[,4]
      q975[1] = bsSip.ci$bca[,5]
      polygon(
        c(min(c(X,-xmax)),min(c(X,-xmax)),
          max(c(X,xmax)),max(c(X,xmax))),
        c(q025[1],q975[1],q975[1],q025[1]),
        col = cols_tr2[4],
        border = NA)
    }
    abline(h = mval[1], col = cols[6], lty=2)
    mtext(text = c('SIP'),
          side = 4,
          at   = mval[1],
          col  = cols[6],
          cex  = 0.75*cex)
  }

  # Delta MUE
  if(showDmue) {
    mval = bsDmue$t0
    if(showCI) {
      q025 = bsDmue.ci$bca[,4]
      q975 = bsDmue.ci$bca[,5]
      polygon(
        c(q025,q975,q975,q025),
        c(0,0,1,1),
        col = cols_tr2[2],
        border = NA)
    }
    abline(v = mval, col=cols[2], lty=2)
    mtext(text = expression(Delta[MUE]),
          side = 3,
          line = 0.6,
          at   = mval,
          col  = cols[2],
          cex  = 0.75*cex)
  }

  box()

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)

}
#' Title
#'
#' @param X
#' @param index
#' @param uX
#' @param ...
#'
#' @return
#'
#' @examples
fsi = function(X, index=1:nrow(X), uX = 0,...){
  # SIP
  v1 = abs(X[index,1])
  v2 = abs(X[index,2])
  N  = length(index)
  if(uX != 0) {
    pert = rnorm(N,0,uX) # Paired datasets
    v1 = v1 + pert
    v2 = v2 + pert
  }
  diff = v1 - v2
  gain = diff < 0 # 1 has smaller errors than 2
  loss = diff > 0
  # How to deal with ties ?
  # eq   = diff == 0
  mg = ifelse(sum(gain) == 0, 0, mean(diff[gain]))
  ml = ifelse(sum(loss) == 0, 0, mean(diff[loss]))
  p  = sum(gain) / N
  # p  = (sum(gain)+0.5*sum(eq)) / N

  return(c(p,mg,ml))
}
#' Title
#'
#' @param X
#' @param index
#' @param uX
#' @param ...
#'
#' @return
#'
#' @examples
dmue = function(X, index=1:nrow(X), uX = 0,...){
  v1 = abs(X[index,1])
  v2 = abs(X[index,2])
  N  = length(index)
  if(uX != 0) {
    pert = rnorm(N,0,uX) # Paired datasets
    v1 = v1 + pert
    v2 = v2 + pert
  }
  return(mean(v1)-mean(v2) )
}