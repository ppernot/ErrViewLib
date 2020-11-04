#' Plot of a set of PDFs
#'
#' @param X
#' @param absErrors
#' @param xlab
#' @param adjust
#' @param xmax
#' @param title
#' @param show.leg
#' @param show.MSE
#' @param col.index
#' @param weights
#' @param units
#' @param label
#' @param leg.lwd
#' @param fill
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotPdf = function(X,
                   absErrors = TRUE,
                   adjust = 1,
                   xlab = NULL,
                   xmin = 0,
                   xmax = NULL,
                   title = '',
                   show.leg = TRUE,
                   show.MSE = FALSE,
                   show.Q95 = TRUE,
                   Q.algo = 'HD',
                   col.index = 1:ncol(X),
                   weights = NULL,
                   units = 'a.u.',
                   label = 0,
                   leg.lwd = 30,
                   fill = FALSE,
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

  if(absErrors)
    X = abs(X)

  if (is.null(xmin))
    xmin = min(X, na.rm = TRUE)

  if (is.null(xmax))
    xmax = max(X, na.rm = TRUE)

  if(is.null(xlab))
    if(absErrors)
      xlab = paste0('|Errors| [',units,']')
  else
    xlab = paste0('Errors [',units,']')

  densx = densy = mu = Q95 = list()
  for (icol in 1:ncol(X)) {
    x = X[, icol]
    sel = !is.na(x)
    x = x[sel]
    mu[[icol]] = mean(x)
    if(Q.algo == 'HD')
      Q95[[icol]] = hd(x, 0.95)
    else
      Q95[[icol]] = quantile(x,0.95,na.rm=TRUE)
    d = density(x, adjust = adjust, kernel = 'rectangular')
    densx[[icol]] = d$x
    densy[[icol]] = d$y
  }

  ylim = c(0, 1.05 * max(unlist(densy)))

  for (icol in 1:ncol(X)) {

    if (icol == 1) {
      plot(
        densx[[icol]], densy[[icol]],
        type = 'l',
        col  = cols[col.index[icol]],
        xlab = xlab,
        xlim = c(xmin, xmax),
        ylim = ylim,
        main = '',
        yaxs = 'i',
        ylab = 'Probability density'
      )
      grid(lwd = 2)
    } else {
      lines(
        densx[[icol]], densy[[icol]],
        type = 'l',
        col  = cols[col.index[icol]]
      )
    }
    if(fill)
      polygon(
        densx[[icol]], densy[[icol]],
        col  = cols_tr[col.index[icol]],
        border = NA
      )

    if(show.Q95)
      abline(v = Q95[[icol]], col = cols[col.index[icol]], lty = 2)

    if (show.MSE)
      abline(v = mu[[icol]], col = cols[col.index[icol]], lty = 3)
  }

  if(!absErrors)
    abline(v=0)

  box()

  if (show.leg & length(colnames(X)) != 0)
    legend(
      'topright',
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