#' Plot Principal Component Analysis results
#'
#' @param X a numeric matrix
#' @param labels a strings vector with the names of the data points
#' @param xlim (optional) a 2-vector of limits for the x-axis
#' @param ylim (optional) a 2-vector of limits for the y-axis
#' @param cex.lab (optional) character expansion factor used for
#'   labelling the points.
#' @param gPars a list of graphical parameters. See \link{setgPars}.
#'
#' @return a plot
#' @export
#'
plotPCA <- function(
  X,
  labels,
  xlim = NULL,
  ylim = NULL,
  cex.lab = 1,
  gPars = ErrViewLib::setgPars()
) {

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  rownames(X)= labels
  pc= prcomp(X, scale.=TRUE)

  par(
    mfrow = c(1, 2),
    mgp = mgp,
    pty = pty,
    mar = c(mar[1:3],2),
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  screeplot(
    pc,
    type ="lines",
    main=' Screeplot',
    col = cols[6]
  )
  grid()
  box()

  if(is.null(xlim) | is.null(ylim)) {
    biplot(
      pc,
      choices=1:2,
      col = c(cols[2],cols[6]),
      cex = cex.lab
    )
  } else {
    # Convert
    biplot(
      pc,
      choices=1:2,
      col = c(cols[2],cols[6]),
      xlim = xlim,
      ylim = ylim,
      cex = cex.lab
    )
  }
  grid()
  box()
}