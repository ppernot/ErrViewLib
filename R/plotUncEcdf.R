#' Plot of a set of ECDFs
#'
#' @param X
#' @param xlab
#' @param xmax
#' @param title
#' @param show.leg
#' @param show.MAE
#' @param col.index
#' @param weights
#' @param units
#' @param label
#' @param leg.lwd
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotUncEcdf = function(X,
                       xlab = NULL,
                       xmax = NULL,
                       title = '',
                       show.leg = TRUE,
                       show.MAE = FALSE,
                       show.Q95 = TRUE,
                       Q.algo = 'HD',
                       col.index = 1:ncol(X),
                       weights = NULL,
                       units = 'a.u.',
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

  if (is.null(xmax))
    xmax = max(X, na.rm = TRUE)

  if(is.null(xlab))
    xlab = paste0('|Errors| [',units,']')

  for (icol in 1:ncol(X)) {
    x = X[, icol]
    sel = !is.na(x)
    x = x[sel]
    io = order(x)
    x = x[io]
    if(!is.null(weights)) {
      prob = cumsum(weights[sel][io])/sum(weights[sel])
    } else {
      prob = (1:length(x)) / length(x)
    }

    if (icol == 1) {
      plot(
        x, prob,
        type = 'l',
        col  = cols[col.index[icol]],
        xlab = xlab,
        xlim = c(0, xmax),
        main = '',
        yaxs = 'i',
        ylab = 'Probability'
      )
      grid(lwd = 2)
    } else {
      lines(x, prob, col = cols[col.index[icol]])
    }

    sigp = sqrt(prob * (1 - prob) / length(x))
    polygon(c(x, rev(x)),
            c(prob - 1.96 * sigp, rev(prob + 1.96 * sigp)),
            col = cols_tr2[col.index[icol]],
            border = NA)

    if(show.Q95) {

      if(Q.algo == 'HD')
        Q95 = hd(x, 0.95)
      else
        Q95 = quantile(x,0.95,na.rm=TRUE)

      segments(Q95, 0, Q95, 0.95,
               col = cols[col.index[icol]], lty = 2)
    }


    if (show.MAE) {
      MAE = mean(abs(x))
      pMAE = prob[which(x >= MAE)[1]]
      segments(MAE, 0, MAE, pMAE,
               col = cols[col.index[icol]], lty = 3)
      segments(0, pMAE, MAE, pMAE,
               col = cols[col.index[icol]], lty = 3)
    }
  }

  if(show.Q95) {
    abline(h = 0.95, col = 2, lty = 2)
    mtext(
      text = '0.95 ',
      at = 0.95,
      side = 2,
      col = 2,
      cex = cex,
      las = 2
    )
  }

  box()

  if (show.leg & length(colnames(X)) != 0)
    legend(
      'bottomright',
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
