#' Generate intervals for local statistics
#'
#' @param X (integer or vector) if integer: size of set to partition,
#'   if vector: ordered set of values to bin
#' @param nBin (integer) number of contiguous intervals
#' @param slide (logical) generate contiguous intervals
#'   (default `slide = FALSE`) or generate intervals from a sliding window
#' @param equiPop (logical) generate intervals  with equal bin counts
#'   (default: `equiPop = TRUE`)
#' @param popMin (integer) minimal bin count in an interval
#' @param logBin (logical) if `equiPop = FALSE`, one can choose between
#'   equal range intervals, or equal log-range intervals
#'   (default `logBin = TRUE`)
#'
#' @return A list with `nbr` the number of intervals, `lwindx` and  `upindx`,
#'   vectors of indices for the first and last point of each interval
#'   `pop` a vector of bin counts, `msg` an error message and the configuration
#'   parameters: `slide`, `equiPop`, `popMin` and `logBin`
#'
#' @export
#'
#' @examples
#' \donttest{
#'   uE  = sqrt(rchisq(1000, df = 4))  # Re-scale uncertainty
#'   int = genIntervals(uE,20) # Generate 20 bins with equal counts
#' }
genIntervals = function(
  X,
  nBin    = 20,
  slide   = FALSE,
  equiPop = TRUE,
  popMin  = 30,
  logBin  = TRUE
) {

  msg = NULL

  if(length(X) == 1) {
    if(is.integer(X) & X > 0)
      N = X
  } else {
    N = length(X)
  }

  if( N / nBin < popMin) {
    # nBin too large vs popMin -> adjust it
    nBin0 = nBin
    nBin  = floor( N / popMin)
    msg   = paste0(
      'nBin was reduced from ', nBin0,' to ',nBin,
      ' in order to respect the minimal bin count,',
      ' popMin = ',popMin
    )
  }

  if (slide) {
    # Sliding interval of width nLoc
    nLoc = floor(N / nBin)
    # Nbr of intervals
    nbr  = N - nLoc +1
    # Lower index of interval in ordered data
    lwindx = 1:nbr
    # Upper index
    upindx = lwindx + nLoc -1
    nBin = nbr

  } else if(equiPop) {
    # Breakpoints for nearly equi-sized intervals
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
    nBin = nbr

    # if(min(upindx-lwindx) < N/nBin/2 |
    #    sum(upindx-lwindx+1) != N      )
    #   stop('>>> Pb in equi-sized intervals design')

  } else {
    # Generate bins with similar ranges and
    # with constraints on min and max populations
    nBin0 = nBin

    xOrd = sort(X)
    xl = xOrd
    if(logBin)
      xl = log(xOrd)
    step = diff(range(xl))/nBin
    lims = min(xl) + (1:nBin -1) * step
    if(logBin)
      lims = exp(lims)

    # Intervals
    lwindx = upindx = c()
    lwindx[1] = 1
    for(i in 2:nBin) {
      lwindx[i] = which(xOrd >= lims[i])[1]
      upindx[i-1] = lwindx[i] - 1
    }
    upindx[nBin] = N

    # Merge small populations
    pop = upindx - lwindx + 1
    while(any(pop < popMin)) {
      # Left tail
      for(i in 1:(nBin-1)) {
        if(pop[i] >= popMin)
          next
        upindx[i]   = upindx[i+1]
        upindx[i+1] = 0
        lwindx[i+1] = 0
        break
      }
      lwindx = lwindx[lwindx!=0]
      upindx = upindx[upindx!=0]
      nBin   = length(lwindx)
      pop    = upindx - lwindx + 1
      # Right tail
      if(pop[nBin] < popMin) {
        upindx[nBin-1] = upindx[nBin]
        lwindx = lwindx[-nBin]
        upindx = upindx[-nBin]
        nBin   = nBin-1
        pop    = upindx - lwindx + 1
      }
    }

    # Split large blobs
    ## Choose upper limit to avoid splitted populations
    ## to become smaller than popMin
    popMax = max(floor(N / nBin0), 2 * popMin)
    while(any(pop > popMax)) {
      for(i in 1:nBin) {
        if(pop[i] <= popMax)
          next
        if(i != nBin) {
          upindx[(i+2):(nBin+1)] = upindx[(i+1):nBin]
          lwindx[(i+2):(nBin+1)] = lwindx[(i+1):nBin]
        }
        upindx[i+1] = upindx[i]
        upindx[i]   = lwindx[i] + floor(pop[i]/2)
        lwindx[i+1] = upindx[i] + 1
        break
      }
      nBin = length(lwindx)
      pop = upindx - lwindx + 1
    }
  }

  intrv = list(
    nbr     = nBin,
    lwindx  = lwindx,
    upindx  = upindx,
    pop     = upindx -lwindx + 1,
    msg     = msg,
    slide   = slide,
    equiPop = equiPop,
    popMin  = popMin,
    logBin  = logBin
  )
  return(intrv)
}