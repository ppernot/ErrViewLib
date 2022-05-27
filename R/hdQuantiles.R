#' Estimate quantiles by Harrell & Davis 1982
#'
#' @param x (vector) a set of values
#' @param q (numeric) a probability level
#' @param na.rm (logical) remove NAs
#'
#' @return A quantile.
#' @export
#'
hd = function (x, q = 0.5, na.rm = TRUE){
  # Extracted from package WRS2 which does not export it.
  # Mair, P., & Wilcox, R. R. (2019). "Robust Statistical Methods
  # in R Using the WRS2", Behavior Research Methods, Forthcoming...
  if (na.rm)
    x = elimna(x)
  n <- length(x)
  m1 <- (n + 1) * q
  m2 <- (n + 1) * (1 - q)
  vec <- seq(along = x)
  w <- pbeta(vec/n, m1, m2) - pbeta((vec - 1)/n, m1, m2)
  y <- sort(x)
  hd <- sum(w * y)
  hd
}
#' Quantiles by Harrell & Davis 1982 for a vector of probabilities
#'
#' @param x (vector) a set of values
#' @param p (vector) a vector of probability levels.
#' @param na.rm (logical) remove NAs.
#'
#' @return A vector of quantiles.
#' @export
#'
vhd = function(x, p=c(0.025,0.975), na.rm = TRUE) {
  # Calculate quantiles by hd for a vector of probabilities
  q = vector(length=length(p))
  for(i in seq_along(p))
    q[i] = hd(x,p[i], na.rm)
  return(q)
}
#' Eliminate NA's
#'
#' @param m
#'
#' @return
#'
elimna = function (m){
  # Used by hd().
  # Extracted from package WRS2 which does not export it.
  # Mair, P., & Wilcox, R. R. (2019). "Robust Statistical Methods
  # in R Using the WRS2", Behavior Research Methods, Forthcoming...
  if (is.null(dim(m)))
    m <- as.matrix(m)
  ikeep <- c(1:nrow(m))
  for (i in 1:nrow(m))
    if (sum(is.na(m[i, ]) >= 1))
      ikeep[i] <- 0
  elimna <- m[ikeep[ikeep >= 1], ]
  elimna
}