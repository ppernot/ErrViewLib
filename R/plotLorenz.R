#' Plot of a set of Lorenz curves
#'
#' @param X
#' @param title
#' @param show.leg
#' @param show.norm
#' @param col.index
#' @param label
#' @param leg.lwd
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotLorenz = function(
  X,
  title = ' Lorenz curves',
  show.norm = FALSE,
  show.leg = TRUE,
  col.index = 1:ncol(X),
  label = 0,
  leg.lwd = 2,
  gPars) {

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

  # gini = lac = c()
  for (icol in 1:ncol(X)) {
    x = X[, icol]
    x = sort(x[ !is.na(x) ])
    prob = (1:length(x)) / length(x)
    lc = cumsum(x)/sum(x)
    # gini[icol] = ineq::Gini(x)
    # lac[icol] = ineq::Lasym(x)

    if (icol == 1) {
      plot(
        prob,lc,
        type = 'l',
        lwd = 2*lwd,
        main = title,
        col  = cols[col.index[icol]],
        xlab = 'p',
        xlim = c(0,1),
        xaxs = 'i',
        ylab = 'L(p)',
        ylim = c(0,1),
        yaxs = 'i'
      )
      grid(lwd = 2)
    } else {
      lines(prob, lc, lwd = 2*lwd, col = cols[col.index[icol]])
    }
  }
  abline(a=0,b=1)
  if(show.norm) {
    x = sort(abs(rnorm(10000,0,1)))
    prob = (1:length(x)) / length(x)
    lc = cumsum(x)/sum(x)
    lines(prob, lc, lwd = 2*lwd, lty=2, col='gray70')
  }
  box()

  if (show.leg & ncol(X) >= 1) {
    # legend = paste0(colnames(X),
    #                 '; Gini = ',signif(gini,2),
    #                 '; LAC = ', signif(lac,2))
    # if(show.norm)
    #   legend = c(legend, 'Normal ; Gini = 0.41; LAC = 0.85')

    legend = colnames(X)
    if(show.norm)
      legend = c(legend, 'Normal')

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
