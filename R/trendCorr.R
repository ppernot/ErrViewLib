#' trendCorr
#'
#' @param Data -
#' @param Errors -
#' @param degree -
#'
#' @return
#' @export
#'
trendCorr <- function(Data, Errors, degree = 0) {

  # Build regression formula
  fo = y ~ 1
  if (degree > 0)
    fo = as.formula(
      paste0('y ~ 1 +',
             paste0(
               'I(x^', 1:degree, ')',
               collapse = '+'
             )))

  # Apply it
  for (i in 1:ncol(Errors)) {
    x = Data[, i]
    y = Errors[, i]
    y = residuals(lm(fo))
    Errors[, i] = y
  }

  # Change names
  prefix = 'c-'
  if(degree == 1)
    prefix = 'lc-'
  else if(degree == 2)
    prefix = 'qc-'
  colnames(Errors) = paste0(prefix, colnames(Errors))

  return(Errors)
}