#' Plot Gini vs. LAC coefs
#'
#' @param X -
#' @param title -
#' @param show.leg -
#' @param show.norm -
#' @param col.index -
#' @param label -
#' @param leg.lwd -
#' @param gPars -
#'
#' @return
#' @export
#'
plotGiniVsLAC = function(
  X,
  title = '',
  show.norm = FALSE,
  show.leg = TRUE,
  col.index = 1:ncol(X),
  label = 0,
  leg.lwd = 2,
  gPars = ErrViewLib::setgPars()
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
    xaxs = 'i'
  )

  if (length(X) == 0)
    return()

  if (class(X) == 'numeric')
    X = as.matrix(X, ncol = 1)

  if (class(X) == 'list') {
    n = names(X)
    X = as.matrix(as.data.frame(X))
    colnames(X) = n
  }

  gini = lac = c()
  for (icol in 1:ncol(X)) {
    x = X[, icol]
    x = sort(x[ !is.na(x) ])
    prob = (1:length(x)) / length(x)
    lc = cumsum(x)/sum(x)
    gini[icol] = ineq::Gini(x)
    lac[icol] = Zanardi(x) #ineq::Lasym(x)
  }

  plot(
    gini, lac,
    pch = 1,
    xlim = c(0,1), #c(0.8*min(c(0.414,gini)),1.2*max(gini)),
    xlab = 'Gini',
    ylim = c(-1,1), #c(0.6,1.4), #c(0.8*min(c(0.85,lac)),1.2*max(lac)),
    ylab = 'Zanardi',
    col = cols[col.index]
  )
  grid(lwd = 2)
  abline(v = 0.414, lty = 2, col = 'gray70') # Ref. Normal
  abline(h = -0.2,  lty = 2, col = 'gray70') # Ref. Normal
  box()

  if (show.leg & ncol(X) >= 1) {
    legend = colnames(X)
    legend(
      'topleft',
      legend = legend,
      bty = 'n',
      col = c(cols[col.index],'gray70'),
      lty = c(rep(1,ncol(X)),2),
      lwd = leg.lwd,
      cex = cex.leg
    )
  }

  if(label >0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)

}
