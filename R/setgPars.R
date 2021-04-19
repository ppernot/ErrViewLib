#' Set default graphical parameters
#'
#' @param type (string) plot type. Default: 'plot'. Other options are 'shiny' and 'publish'.
#' @param gPars (list) list of graphical parameters to supersede defaults
#'
#' @return A list with the following arguments:
#' \describe{
#'   \item{cols}{A set of colors}
#'   \item{cols_tr}{A set of highly transparent colors}
#'   \item{cols_tr2}{Sames as cols_tr with more density}
#'   \item{pty}{A character specifying the type of plot region to be used;
#'     's' generates a square plotting region and 'm' generates the maximal
#'     plotting region (default: 's')}
#'   \item{mar}{Margins around the plot. The default is c(3,3,3,.5).
#'     See \link{graphics::par}.
#'   \item{mgp}{The margin line (in mex units) for the axis title,
#'     axis labels and axis line. The default is c(2,.75,0).}
#'   \item{tcl}{The length of tick marks as a fraction of the height
#'     of a line of text. The default is -0.5.}
#'   \item{lwd}{Width of the plotted lines. The default is 4.}
#'   \item{cex}{A numerical value giving the amount by which
#'     plotting text and symbols should be magnified relative
#'     to the default. The default value is 4.}
#'   \item{cex.leg}{A numerical value giving the amount by which
#'     the legends should be magnified relative
#'     to the default. The default value is 0.7.}
#'   \item{reso}{A numerical value giving the nominal pixel size
#'     for a PNG plot. The default is 1200.}
#' }
#'
#' @export
#'
#' @examples
setgPars = function(
  type = c('plot','publish','shiny'),
  gPars = list()
) {

  type = match.arg(type)

  # Defaults
  gParsDef = list(
    cols     = rev(inlmisc::GetColors(8))[1:7],
    cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.2))[1:7],
    cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.5))[1:7],
    pty      = 's',
    mar      = c(3,3,3,.5),
    mgp      = c(2,.75,0),
    tcl      = -0.5,
    reso     = 1200
  )
  # Adapt to plot type
  if(type == 'plot') {
    gParsDef$lwd     = 1.5
    gParsDef$cex     = 1
    gParsDef$cex.leg = 0.7

  } else if(type == 'shiny') {
    gParsDef$mar     = c(3,3,1,1)
    gParsDef$lwd     = 2
    gParsDef$cex     = 1.4
    gParsDef$cex.leg = 0.8

  } else {
    # 'publish' (for png plots with base resolution 'reso')
    gParsDef$lwd     = 4
    gParsDef$cex     = 4
    gParsDef$cex.leg = 0.7
  }

  # Override by user's specs, if any
  if(!is.list(gPars)) {
    message('>>> Argument gPars should be a list !\n')
    stop(call.=FALSE)
  } else {
    if(length(gPars) != 0) {
      for (n in names(gPars)) {
        if(!n %in% names(gParsDef)) {
          message('>>> Incorrect variable name in gPars: ',n,'\n')
          stop(call.=FALSE)
        }
        gParsDef[[n]] = rlist::list.extract(gPars, n)
      }
    }
  }

  return(gParsDef)
}


