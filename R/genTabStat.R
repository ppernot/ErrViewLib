#' Generate results table from benchmarking statistics generated
#' by 'estBS()'
#'
#' @param S
#' @param comp
#' @param ref
#' @param numDig (integer) number of digits to keep for uncertainty display.
#' @param units (string) units of the data.
#' @param short (logical) use parenthetic notation to display uncertainty.
#'
#' @return A dataframe.
#' @export
#'
#' @examples
genTabStat = function(
  S,
  comp    = TRUE,
  ref     = 0,
  numDig  = 1,
  units   = 'a.u.',
  short   = TRUE
) {

  colUnc = function(
    prop,
    x,
    ux,
    units = '',
    short = TRUE,
    numDig=2
  ) {
    # Generate matrix of values and uncertainties with adequate truncation
    # and format
    if(short) {
      vu = matrix(
        apply(cbind(x,ux),1,
              function(x) prettyUnc(x[1],x[2],numDig = numDig)),
        ncol=1)
      vu = rbind(units,vu)
      colnames(vu) = prop
    } else{
      vu = matrix(
        apply(cbind(x,ux),1,
              function(x) unlist(formatUnc(x[1],x[2],numDig = numDig))),
        ncol=2,
        byrow = TRUE)
      vu = rbind(c(units,units),vu)
      colnames(vu) = c(prop, paste0('u_', prop))
    }
    return(vu)
  }
  methods = names(S[[S[['props']][1]]]$val)
  nm = length(methods)

  df = data.frame(Methods = c('Units',methods)) # Leave 1rst row for units
  for (prop in S[['props']]) {

    un = units
    if(prop %in% c('P1','W',
                   'skew','kurt','skewgm','kurtcs',
                   'gini','gimc','pietra','lasym','Zanardi')
    )
      un = ''

    v  = S[[prop]]$val
    uv = S[[prop]]$unc
    vu = colUnc(prop,v,uv,units=un,short=short,numDig=numDig)
    df = cbind(df, vu)

    if (comp &
        (prop %in% c('mue','rmsd','q95hd')) &
        nm > 1) {
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
        udiff = sqrt(S[[prop]]$unc[mi]^2 + S[[prop]]$unc[mj]^2)
        compt[j] = 2*(1-pnorm(diff/udiff))
      }
      compt = matrix(compt,ncol=1)
      colnames(compt) = paste0('punc_',prop)
      df = cbind(df,rbind('',compt))

      # t-test for paired values
      compt = c()
      compt[im] = 1
      for (j in (1:nm)[-im]) {
        mj = methods[j]
        diff  = unlist(S[[prop]]$bs[mi]) - unlist(S[[prop]]$bs[mj])
        compt[j] = genpval(diff)
      }
      compt = matrix(compt,ncol=1)
      colnames(compt) = paste0('pg_',prop)
      df = cbind(df,rbind('',compt))

      # Pinv
      compt = rep(NA,nm)
      for (j in (1:nm)[-im]) {
        mj    = methods[j]
        d0    = S[[prop]]$val[mi] - S[[prop]]$val[mj]
        diff  = unlist(S[[prop]]$bs[mi]) - unlist(S[[prop]]$bs[mj])
        compt[j] = round(pinv(diff,d0),2)
      }
      compt = matrix(compt,ncol=1)
      colnames(compt) = paste0('Pinv_',prop)
      df = cbind(df,rbind('',compt))
    }
  }

  if(!is.null(S$sip)) {
    # Mean SIP
    prop = 'MSIP'
    msip = rowMeans(S[['sip']], na.rm = TRUE)
    umsip = sqrt(rowSums(S[['usip']] ^ 2, na.rm = TRUE) / nm) # Hyp. indép.
    vu = colUnc(prop, msip, umsip, units = '', short = short, numDig = numDig)
    df = cbind(df, vu)

    # SIP for best MUE
    prop = 'SIP'
    mi   = methods[which.min(S[['mue']]$val)]
    sip  = S[['sip']][mi, ]
    usip = S[['usip']][mi, ]
    vu = colUnc(prop, sip, usip, units = '', short = short, numDig = numDig)
    df = cbind(df, vu)

    # Mean gain
    prop = 'MG'
    mg  = S[['mg']][mi,]
    umg = S[['umg']][mi,]
    vu = colUnc(prop, mg, umg, units = units, short = short, numDig = numDig)
    df = cbind(df, vu)

    # Mean loss
    prop = 'ML'
    mg = -S[['mg']][, mi]
    umg = S[['umg']][, mi]
    vu = colUnc(prop, mg, umg, units = units, short = short, numDig = numDig)
    df = cbind(df, vu)
  }
  return(df)
}
#' Generalized p-value (Liu & Sing 1997, Wilcox 2012)
#'
#' @param X (vector) a set of values
#'
#' @return a p-value
#' @export
#'
#' @examples
genpval = function(X) {
  ps = (sum(X < 0) + 0.5 * sum(X == 0)) / length(X)
  2 * min(ps, 1 - ps)
}
#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
pval = function(x) {
  2*pnorm(x, lower.tail = FALSE)
}
#' Probability to have a sign different from the sign of d0
#'
#' @param X (vector) values to be tested.
#' @param d0 (numeric) reference value.
#'
#' @return A probability.
#' @export
#'
#' @examples
pinv = function (X,d0) {
  # The zeros (sign = 0) are subtracted
  A = sum( sign(X) != sign(d0) )
  C = sum( X == 0 )
  (A - C)/length(X)
}
#' Truncate value and uncertainty to consistent number of digits.
#'
#' @param y (numeric) value
#' @param uy (numeric) uncertainty on `y`
#' @param numDig (numeric) number of digits to keep on `uy`
#'
#' @return A list with strings of truncated values of `y` and `uy`.
#' @export
#'
#' @examples
formatUnc = function(y, uy, numDig = 2) {

  if (!is.finite(y) | !is.finite(uy) | uy <= 0)
    return(
      list(y  = y, uy = uy)
    )

  # Get scales
  n0 = 1 + floor(log10(abs(y)))
  n1 = floor(log10(uy))

  # Format uncertainty
  fmt = switch(
    sign(n1) + 2, # Map (-1,0,1) to (1,2,3)
    paste0("%", n0 - n1 + numDig - 1, ".", -n1 + numDig - 1, "f"),
    paste0("%", n0 - n1 + numDig - 1, ".", -n1 + numDig - 1, "f"),
    paste0("%", n0, ".0f")
  )
  short_y  = sprintf(fmt, y)
  short_uy = sprintf(fmt,uy) #paste0(signif(uy, numDig))

  return(
    list(
      y  = short_y,
      uy = short_uy
    )
  )
}
#' Print value and uncertainty in parenthesis format
#'
#' @param y (numeric) value
#' @param uy (numeric) uncertainty on `y`
#' @param numDig (numeric) number of digits to keep on `uy`
##'
#' @return A string.
#' @export
#'
#' @examples
prettyUnc = function(y, uy, numDig = 2) {

  if (!is.finite(uy) |
      !is.finite(y)  |
      is.na(y)       |
      is.na(uy)      |
      uy <0            )
    return(NA)

  # Get scales
  n0 = 1 + floor(log10(abs(y)))
  n1 = floor(log10(uy))

  # Format uncertainty
  switch(
    sign(n1) + 2, # Map (-1,0,1) to (1,2,3)
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
    }
  )
  short_y  = sprintf(fmt, y)

  str   = paste0(short_y, '(', short_uy, ')')

  return(str)
}
