#' Title
#'
#' @param S
#' @param comp
#' @param ref
#' @param numDig
#'
#' @return
#' @export
#'
#' @examples
genTabStat = function(S, comp=TRUE, ref = 0, numDig=1) {
  methods = names(S[[S[['props']][1]]]$val)
  nm = length(methods)

  df = data.frame(Methods = methods)
  for (prop in S[['props']]) {
    v = S[[prop]]$val
    vu = matrix(
      apply(cbind(v,S[[prop]]$unc),1,
            function(x) prettyUnc(x[1],x[2],numDig = numDig)),
      ncol=1)
    colnames(vu) = prop
    df = cbind(df, vu)
    if (comp & (prop %in% c('mue','wmue','rmsd','q95hd')) ) {
      if(ref != 0){
        # compare with specified method
        im = ref
      } else {
        # compare with best score
        im = which.min(abs(v))
      }
      mi = methods[im]

      # t-test for unpaired  values
      compt = c()
      compt[im] = 1
      for (j in (1:nm)[-im]) {
        mj = methods[j]
        diff  = abs(S[[prop]]$val[mi] - S[[prop]]$val[mj])
        udiff = sqrt(
          S[[prop]]$unc[mi]^2 + S[[prop]]$unc[mj]^2
        )
        compt[j] = 2*(1-pnorm(diff/udiff))
      }
      names(compt) = methods
      df = cbind(df, punc = round(compt,3))

      # t-test for paired values
      compt = c()
      compt[im] = 1
      for (j in (1:nm)[-im]) {
        mj = methods[j]
        diff  = unlist(S[[prop]]$bs[mi]) -
          unlist(S[[prop]]$bs[mj])
        compt[j] = genpval(diff)
      }
      names(compt) = methods
      df = cbind(df, pg = round(compt,3))

      # Pinv
      compt = c()
      compt[im] = NA
      for (j in (1:nm)[-im]) {
        mj = methods[j]
        d0    = S[[prop]]$val[mi] - S[[prop]]$val[mj]
        diff  = unlist(S[[prop]]$bs[mi]) -
          unlist(S[[prop]]$bs[mj])
        compt[j] = pinv(diff,d0)
      }
      names(compt) = methods
      df = cbind(df, Pinv = round(compt,3))
    }
  }

  if(!is.null(S$sip)) {
    # Mean SIP
    msip = rowMeans(S[['sip']], na.rm=TRUE)
    umsip = sqrt(rowSums(S[['usip']]^2, na.rm=TRUE) / nm) # Hyp. indép.
    vu = matrix(
      apply(cbind(msip,umsip),1,
            function(x) prettyUnc(x[1],x[2],numDig = numDig)),
      ncol=1)
    colnames(vu) = 'MSIP'
    rownames(vu) = methods
    df = cbind(df, vu)

    # SIP for best MUE
    v = S[['mue']]$val
    im = which.min(abs(v))
    mi = methods[im]
    sip = S[['sip']][mi,]
    usip = S[['usip']][mi,]
    vu = matrix(
      apply(cbind(sip,usip),1,
            function(x) prettyUnc(x[1],x[2],numDig = numDig)),
      ncol=1)
    colnames(vu) = 'SIP'
    rownames(vu) = methods
    df = cbind(df, vu)

    # Mean gain
    mg  = S[['mg']][mi,]
    umg = S[['umg']][mi,]
    vu = matrix(
      apply(cbind(mg,umg),1,
            function(x) prettyUnc(x[1],x[2],numDig = numDig)),
      ncol=1)
    colnames(vu) = 'MG'
    rownames(vu) = methods
    df = cbind(df, vu)

    # Mean loss
    mg = -S[['mg']][,mi]
    umg = S[['umg']][,mi]
    vu = matrix(
      apply(cbind(mg,umg),1,
            function(x) prettyUnc(x[1],x[2],numDig = numDig)),
      ncol=1)
    colnames(vu) = 'ML'
    rownames(vu) = methods
    df = cbind(df, vu)

  }

  return(df)
}
#' Title
#'
#' @param X
#'
#' @return
#'
#' @examples
genpval = function(X) {
  # Generalized p-value (Liu & Sing 1997, Wilcox 2012)
  ps = (sum(X < 0) + 0.5 * sum(X == 0)) / length(X)
  2 * min(ps, 1 - ps)
}
#' Title
#'
#' @param X
#' @param d0
#'
#' @return
#'
#' @examples
pinv = function (X,d0) {
  # Probability to have a sign different from d0's
  # The zeros (sign = 0) are excluded
  A = sum( sign(X) != sign(d0) )
  C = sum( X == 0 )
  (A - C)/length(X)
}
#' Title
#'
#' @param y
#' @param uy
#' @param numDig
#'
#' @return
#'
#' @examples
prettyUnc = function(y, uy, numDig = 2) {
  # Print result + uncertainty in parenthesis format

  if (!is.finite(y))
    return(y)

  if (!is.finite(uy) | uy<=0)
    return(y)

  # Get scales
  n0 = 1 + floor(log10(abs(y)))
  n1 = floor(log10(uy))

  # Format uncertainty
  switch(sign(n1) + 2, # Map (-1,0,1) to (1,2,3)
         {
           fmt = paste0("%", n0 - n1 + numDig - 1, ".", -n1 + numDig - 1, "f")
           short_uy = signif(uy / 10 ^ (n1 - numDig + 1), numDig)
         },
         {
           fmt = paste0("%", n0 - n1 + numDig - 1, ".", -n1 + numDig - 1, "f")
           short_uy = signif(uy / 10 ^ n1, numDig)
         },
         {
           fmt = paste0("%", n0, ".0f")
           short_uy = signif(uy / 10 ^ (n1 - numDig + 1), numDig)
         })
  short_y  = sprintf(fmt, y)

  return(paste0(short_y, '(', short_uy, ')'))
}