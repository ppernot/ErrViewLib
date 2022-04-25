#' Plot of Systematic Improvement Probability matrix
#'
#' Interface to \code{corrplot::corrplot()}
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
plotSIPMat = function(
  X,
  method = 'circle',
  order = TRUE,
  main = '',
  label  = 0,
  cex.lab = 1,
  gPars = ErrViewLib::setgPars()
){

  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    # pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  col2 <- colorRampPalette(
    c("#67001F", "#B2182B", "#D6604D", "#F4A582",
      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
      "#4393C3", "#2166AC", "#053061"))
  colors = rev(c(col2(200),col2(200))) # Bypass bug in corrplot

  diag(X) = 0 # Replace NAs
  if (order) {
    io = order(rowMeans(X),decreasing = TRUE)
    X = X[io,io]
  }
  corrplot::corrplot(
    X,
    cl.lim = c(0,1),
    col = colors,
    is.corr = FALSE,
    diag = TRUE,
    method = method,
    order = 'original',
    tl.col = 'black',
    tl.cex = cex.lab,
    cl.cex = cex.lab)

  title(main,line=0)

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj  = 0.975,
      # las  = 1,
      cex  = cex,
      line = 0.2)

  # if(unc) {
  #   ncolors=11
  #   colors=fields::two.colors(
  #     n=ncolors, start="white", end="red",middle="pink")
  # } else {
  #   ncolors = 11
  #   colors=fields::two.colors(
  #     n=ncolors, start="blue", end="red", middle="white")
  # }
  # # colors = inlmisc::GetColors(11,scheme ='iridescent')
  # X[is.na(X)] = -0.05
  # zlim = c(-0.1,1)
  #

}
