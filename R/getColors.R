#' Generate colors according to smooth rainbow scheme
#'
#' @param n (integer) number of colors
#' @param alpha (numeric) level of alpha channel in [0,1]
#'
#' @return a color palette
#'
myGetColors = function (n,alpha = NULL) {
  color <- c("#E8ECFB","#DDD8EF","#D1C1E1","#C3A8D1","#B58FC2","#A778B4","#9B62A7","#8C4E99","#6F4C9B","#6059A9","#5568B8","#4E79C5","#4D8AC6","#4E96BC","#549EB3","#59A5A9","#60AB9E","#69B190","#77B77D","#8CBC68","#A6BE54","#BEBC48","#D1B541","#DDAA3C","#E49C39","#E78C35","#E67932","#E4632D","#DF4828","#DA2222","#B8221E","#95211B","#721E17","#521A13")
  names(color) <- paste0('color',1:length(color))
  value = seq_along(color)
  value <- scales::rescale(value)
  x <- seq.int(0, 1, length.out = 255)
  color <- (scales::gradient_n_pal(color, values = value))(x)
  pal <- (grDevices::colorRampPalette(color, bias = 1, space = "Lab"))(n)
  if (!is.null(alpha))
    pal <- grDevices::adjustcolor(pal, alpha.f = alpha)
  return(pal)
}