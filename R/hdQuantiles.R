#' Title
#'
#' @param x
#' @param q
#' @param na.rm
#'
#' @return
#' @export
#'
#' @examples
hd = function (x, q = 0.5, na.rm = TRUE){
  # Estimate quantiles by Harrell & Davis 1982
  # Extracted from package WRS2 which does not export it
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
#' Title
#'
#' @param m
#'
#' @return
#'
#' @examples
elimna = function (m){
  if (is.null(dim(m)))
    m <- as.matrix(m)
  ikeep <- c(1:nrow(m))
  for (i in 1:nrow(m))
    if (sum(is.na(m[i, ]) >= 1))
      ikeep[i] <- 0
  elimna <- m[ikeep[ikeep >= 1], ]
  elimna
}