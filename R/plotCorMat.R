#' Title
#'
#' @param X
#' @param method
#' @param order
#' @param main
#' @param label
#' @param cex.lab
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotCorMat = function(X,
                      method = "ellipse",
                      order  = "original",
                      main  = '',
                      label  = 0,
                      cex.lab = 1,
                      gPars){

  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  M = corrplot::corrplot(
    X,
    method = method,
    order = order,
    tl.col = 'black',
    tl.cex = cex.lab)

  mtext(main, line=0.5, cex = 1.25*cex, adj = 0.025)

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj  = 0.975,
      # las  = 1,
      cex  = cex,
      line = 0.15)

  orderedNames = colnames(M)
  invisible(orderedNames)
}
