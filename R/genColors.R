#' Color palette for parallel plots
#'
#' @param sample
#'
#' @return
#'
genColors = function(sample) {
  ncols=length(sample)
  co=(    sample -min(sample))/
    (max(sample)-min(sample))
  indx=round(1+(ncols-1)*co)
  # cols=fields::two.colors(ncols,start="blue",middle="yellow",end="red")[indx]
  cols = inlmisc:: GetColors(ncols,scheme = 'sunset', alpha=0.8)[indx]
  return(cols)
}
#' Color palette for parallel plots
#'
#' @param sample
#'
#' @return
#'
genColorsSample = function(sample) {
  ncols=128
  co=(    sample -min(sample))/
    (max(sample)-min(sample))
  indx=round(1+(ncols-1)*co)
  cols=fields::two.colors(ncols,start="gold",middle="white",end="purple")[indx]
  return(cols)
}
