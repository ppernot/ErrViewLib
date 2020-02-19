#' N-out of-N Bootstrap of order statistics
#'
#' @param E
#' @param score
#' @param nMC
#'
#' @return
#' @export
#'
#' @examples
rankBS = function(E, score = 'mue', nMC = 1000) {
  # Bootstrap rank statistics
  # (based on boot::boot)

  if(score == 'msip') {
    bs = boot::boot(
      E,
      statistic = fRankMSIP,
      R = nMC
    )
  } else {
    bs = boot::boot(
      E,
      statistic = frank,
      R = nMC,
      fscore = get(score)
    )
  }

  mRank = bs$t0
  pRank = round(
    apply( bs$t, 2,
           function(x) {
             f = rep(0, ncol(E))
             for(i in 1:length(f))
               f[i] = mean(x == i)
             return(f)
           }
    ),
    2
  )
  rownames(pRank) = colnames(E)
  colnames(pRank) = paste0(1:ncol(E))

  return(list(mRank = mRank, pRank=pRank))

}
#' M-out of-N Bootstrap of order statistics
#'
#' @param E
#' @param score
#' @param nMC
#' @param M
#'
#' @return
#' @export
#'
#' @examples
rankBS2 = function(E, score = 'mue', nMC = 1000, M = nrow(E)) {
  # Bootstrap rank statistics
  # Option for M-outof-N bs (based on distillery::booter)

  if(score == 'msip') {
    bs = distillery::booter(
      x = E,
      statistic = fRankMSIP,
      B = nMC,
      rsize = M
    )
  } else {
    bs = distillery::booter(
      x = E,
      statistic = frank,
      B = nMC,
      rsize = M,
      fscore = get(score)
    )
  }

  mRank = bs$original.est
  pRank = round(
    apply( t(bs$results), 2,
           function(x) {
             f = rep(0, ncol(E))
             for(i in 1:length(f))
               f[i] = mean(x == i)
             return(f)
           }
    ),
    2
  )
  rownames(pRank) = colnames(E)
  colnames(pRank) = paste0(1:ncol(E))

  return(list(mRank = mRank, pRank=pRank))

}
#' Results table for ranking probabilities
#'
#' @param E
#' @param scores
#' @param nMC
#'
#' @return
#'
#' @examples
tabRank = function (E, scores=c('mue','q95hd'), nMC=1000) {
  r = list()
  for (score in scores) {
    tab = rankBS(E,score,nMC)
    mrank = order(tab$mRank)
    prank = apply(tab$pRank,1,max)
    rank  = apply(tab$pRank,1,which.max)
    r[[score]] = data.frame(mrank,rank,prank)
  }
  return(r)
}
