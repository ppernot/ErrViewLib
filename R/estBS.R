#' Estimate benchmarking statistics and their uncertainty by bootstrap
#'
#' @param error
#' @param props
#' @param nboot
#' @param do.sip
#' @param eps
#' @param seed
#' @param silent
#'
#' @return
#' @export
#'
#' @examples
estBS1 = function(error,
                  props = c('mue',
                            'mse',
                            'rmsd',
                            'skew',
                            'kurt',
                            'q95',
                            'q95hd',
                            'P1',
                            'W'),
                  nboot  = 1000,
                  do.sip = TRUE,  # Generate SIP stats
                  eps    = 1,     # Threshold for P1
                  seed   = 123,
                  silent = TRUE) {

  # Process data
  results = list()
  methods = names(error)
  nm = length(methods)
  results[['props']] = props
  for (prop in props) {
    results[[prop]] = list()
    if(!silent)
      print(prop)

    statistic = get(prop) # associated function
    bsl = list()
    for (i in 1:nm) {
      m = methods[i]
      set.seed(seed) # Correlate bs for correl. estim.
      bs = boot::boot(error[[m]], statistic, R=nboot, eps = eps)
      bsl[[m]] = bs
      results[[prop]][['val']][[m]] = bs$t0
      results[[prop]][['unc']][[m]] = sd(bs$t)
      results[[prop]][['bs' ]][[m]] = bs$t
    }

    if(nm > 1) {
      # Correlation of scores
      C = matrix(1, nrow = nm, ncol = nm)
      colnames(C) = rownames(C) = methods
      for (i in 1:(nm - 1)) {
        mi = methods[i]
        for (j in (i + 1):nm) {
          mj = methods[j]
          C[i, j] = cor(bsl[[mi]]$t, bsl[[mj]]$t)
          C[j, i] = C[i, j]
        }
      }
      results[[prop]][['corr']] = C
    }
  }

  # Systematic improvement probability
  if(do.sip & nm > 1) {
    fsi = function(X, index=1:nrow(X), uX = 0,...){
      v1 = abs(X[index,1])
      v2 = abs(X[index,2])
      N  = length(index)
      if(uX != 0) {
        pert = rnorm(N,0,uX) # Paired datasets
        v1 = v1 + pert
        v2 = v2 + pert
      }
      diff = v1 - v2
      gain = diff < 0 # 1 has smaller errors than 2
      loss = diff > 0
      # How to deal with ties ?
      # eq   = diff == 0
      mg = ifelse(sum(gain) == 0, 0, mean(diff[gain]))
      # ml = ifelse(sum(loss) == 0, 0, mean(diff[loss]))
      p  = sum(gain) / N
      # p  = (sum(gain)+0.5*sum(eq)) / N

      return(c(p,mg))
    }
    sip = usip = mg = umg = matrix(NA, nrow = nm, ncol = nm)
    for (i in 1:nm) {
      mi = methods[i]
      for (j in 1:nm) {
        if (j==i) next
        mj = methods[j]
        bs = boot::boot(
          cbind(error[[mi]],error[[mj]]),
          statistic = fsi,
          R=nboot,
          uX = 0)
        sip[i, j] = bs$t0[1]
        usip[i,j] = sd(bs$t[,1],na.rm=TRUE)
        mg[i, j]  = bs$t0[2]
        umg[i,j]  = sd(bs$t[,2],na.rm=TRUE)
      }
    }
    rownames(sip) =
      rownames(usip) =
      rownames(mg) =
      rownames(umg) =
      colnames(sip) =
      colnames(usip) =
      colnames(mg) =
      colnames(umg) = methods
    results[['sip']]  = sip
    results[['usip']] = usip
    results[['mg']]   = mg
    results[['umg']]  = umg
  }

  return(results)
}
