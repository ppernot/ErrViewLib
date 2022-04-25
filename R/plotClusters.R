#' Plot dendrograms for hierarchical clusering of systems or methods
#'
#' @param X a numeric matrix of error values
#' @param type (optional) a string describing if "systems" or "methods"
#'    are to be clustered. Default : "systems".
#' @param method (optional) a string describing the hierarchical
#'    clustering method. Default : "complete".
#' @param xlim (optional) a 2-vector of limits for the x-axis.
#' @param ylim (optional) a 2-vector of limits for the y-axis.
#' @param cex.lab (optional) a scale factor for the size of labels.
#' @param gPars a list of graphical parameters. See \link{setgPars}
#'
#' @return a plot
#' @export
#'
#' @examples
plotClusters <- function(
  X,
  type   = c("systems","methods"),
  method = c("complete","single"),
  xlim = NULL,
  ylim = NULL,
  cex.lab = 1,
  gPars = ErrViewLib::setgPars()
) {

  type   = match.arg(type)
  method = match.arg(method)

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  if(type == "methods")
    X = t(X)

  cl = hclust(dist(X), method = method)
  hcd = as.dendrogram(cl)

  par(
    mfrow = c(1,1),
    mgp = mgp,
    pty = "m",
    mar = mar,
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  plot(
    hcd,
    type  = "rectangle",
    hang  = -1,
    check = FALSE,
    nodePar = list(pch = NA, lab.cex = cex.lab, col=cols[1]),
    edgePar = list(col = cols[1], t.col=cols[1]),
    xlab  = toTitleCase(type),
    ylab  = "Height",
    xlim  = xlim,
    ylim  = ylim,
    main  = 'Hierarchical clustering dendrogram'
  )
  rect.hclust(
    cl,
    k = 2,
    border = cols[-1]
  )

}