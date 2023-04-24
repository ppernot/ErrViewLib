#' Tests the fraction of intervals in agreement with target statistics
#' following LCP or LZV/LZISD analysis
#'
#' @param stats An object returned by a LCP or LZV/LZISD analysis
#' @param prob (numeric) the probability levels of the intervals to be tested
#' @param silent (logical) should the function run silently ?
#'
#' @return Invisibly returns a list with binomial test results
#'
#' @export
#'
testCovIntProb <- function(
  stats,
  prob = 0.95,
  silent = FALSE
) {

  if(!is.null(stats$prob)) {
    # LCP analysis
    for(i in seq_along(stats$prob)) {
      p = stats$prob[i]
      fail = 0
      Ns = ncol(stats$pc)
      for(j in 1:Ns) {
        if( (p-stats$pcl[i,j]) *
            (p-stats$pcu[i,j]) > 0)
          fail = fail + 1
      }
      test = binom.test(Ns-fail,Ns, prob)
    }

    invisible(
      list(

      )
    )

  } else {
    # LZV/LZISD analysis
    p = 1
    fail = 0
    Ns = length(stats$pc)
    for(j in 1:Ns) {
      if( (p-stats$pcl[j]) *
          (p-stats$pcu[j]) > 0)
        fail = fail + 1
    }
    test = binom.test(Ns-fail,Ns, prob)
    msg = paste('Valid Z freq.:',round(test$estimate,2),
                paste0('[',round(test$conf.int[1],2),'-',round(test$conf.int[2],2),'];'),
                'Target ~',round(prob,2))
    if(!silent)
      cat(msg,'\n')

    invisible(
      list(
        success  = test$statistic,
        trials   = test$parameter,
        estimate = test$estimate,
        conf.int = test$conf.int,
        p.value  = test$p.value,
        prob     = prob,
        msg      = msg
      )
    )

  }

}