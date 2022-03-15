#' Generate intervals for local statistics
#'
#' @param N (integer) size of set to partition
#' @param nBin (integer) number of contiguous intervals
#' @param slide (logical) generate contiguous intervals
#'   (default `slide = FALSE`) or gerenate from a sliding window
#'
#' @return A list with `nbr` the number of intervals, `lwindx` and  `upindx`,
#'   vectors of indices for the first and last point of each interval.
#'
#' @export
#'
#' @examples
genIntervals = function(
  N,
  nBin,
  slide = FALSE
) {
  if (slide) {
    # Sliding interval of width nLoc
    nLoc = floor(N / nBin)
    # Nbr of intervals
    nbr  = N - nLoc +1
    # Lower index of interval in ordered data
    lwindx = 1:nbr
    # Upper index
    upindx = lwindx + nLoc -1

  } else {
    # Breakpoints of nearly equi-sized intervals
    X    = 1:N
    p    = seq(0, 1, length.out = nBin + 1)[1:nBin]
    br   = ErrViewLib::vhd(X, p = p)
    # Nbr of intervals
    nbr  = length(br)
    # Lower index of interval in ordered data
    lwindx = upindx = c()
    lwindx[1] = 1
    for (i in 2:nbr)
      lwindx[i] = which(X > br[i])[1]
    # Upper index
    for (i in 1:(nbr-1))
      upindx[i] = lwindx[i+1]-1
    upindx[nbr] = N

    if(min(upindx-lwindx) < N/nBin/2 |
       sum(upindx-lwindx+1) != N      )
      stop('>>> Pb in equi-sized intervals design')
  }

  return(
    list(
      nbr = nbr,
      lwindx = lwindx,
      upindx = upindx
    )
  )
}