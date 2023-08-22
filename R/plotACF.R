#' Wrapper for the acf function
#'
#' @param X (vector) a numerical vector
#' @param lag.max (integer)
#' @param title (string)
#' @param col (integer) color index of the CI
#' @param label (integer) index of letter for subplot tag
#' @param gPars (list) graphical parameters
#'
#' @return a plot
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   E   = rnorm(uE, mean=0, sd=uE)    # Generate errors
#'   res = plotLZMS(uE, E/uE, ylim = c(0,2))
#'   plotACF(res$pc)
#' }
plotACF = function(
  X,
  lag.max   = NULL,
  title     = '',
  col       = 5,
  label     = 0,
  gPars     = ErrViewLib::setgPars()
){
  if(length(gPars) == 0)
    gPars = ErrViewLib::setgPars()
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))
  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = 's',
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    cex.main = 1
  )

  acf(
    X,
    lag.max = lag.max,
    main = title,
    ci.col = cols[col]
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