#' Plot of a set of Lorenz curves
#'
#' @param X
#' @param title
#' @param show.leg
#' @param col.index
#' @param label
#' @param leg.lwd
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotUncEcdf = function(X,
                       title = 'Lorenz curves',
                       show.leg = TRUE,
                       col.index = 1:ncol(X),
                       label = 0,
                       leg.lwd = 30,
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


  for (icol in 1:ncol(X)) {
    x = X[, icol]
    sel = !is.na(x)
    x = x[sel]
    io = order(x)
    x = x[io]
    prob = (1:length(x)) / length(x)
    lc = cumsum(x)/sum(x)

    if (icol == 1) {
      plot(
        prob,lc,
        type = 'l',
        col  = cols[col.index[icol]],
        xlab = 'p',
        xlim = c(0,1),
        main = '',
        yaxs = 'i',
        ylab = 'L'
      )
      grid(lwd = 2)
    } else {
      lines(prob, lc, col = cols[col.index[icol]])
    }
    abline(a=0,b=1)

  }
  box()

  if (show.leg & length(colnames(X)) != 0)
    legend(
      'topleft',
      title = title,
      title.adj = 0,
      legend = colnames(X),
      bty = 'n',
      col = cols_tr2[col.index],
      lty = 1,
      lwd = leg.lwd,
      cex = cex.leg
    )

  if(label >0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)

}
