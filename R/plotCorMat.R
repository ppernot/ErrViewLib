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
                      main  = NULL,
                      label  = 0,
                      cex.lab = 1,
                      gPars){

  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(
    mfrow = c(1, 1),
    # mar = mar,
    mgp = mgp,
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  if(is.null(main))
    mar = c(0,0,0,0)

  M = corrplot::corrplot(
    X,
    method = method,
    order = order,
    tl.col = 'black',
    tl.cex = cex.lab,
    cl.cex = 0.75*cex.lab,
    mar = mar)

  if(!is.null(main))
    mtext(main, line=2.7,
          cex = 1.5*cex.lab*cex,
          font = 2,
          adj = 0.2)

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj  = 0.975,
      # las  = 1,
      cex  = 1.5*cex.lab*cex,
      line = ifelse(is.null(main),0.5,2.7)
      )

  orderedNames = colnames(M)
  invisible(orderedNames)
}
