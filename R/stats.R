#' Title
#'
#' @param X
#' @param index
#'
#' @return
#'
#' @examples
mse = function(X, index = 1:length(X)) {
  mean(X[index])
}
#' Title
#'
#' @param X
#' @param index
#'
#' @return
#'
#' @examples
rmsd = function(X, index = 1:length(X)) {
  sd(X[index])
}
#' Title
#'
#' @param X
#' @param index
#'
#' @return
#'
#' @examples
mue = function(X, index = 1:length(X)) {
  mean(abs(X[index]))
}
#' Title
#'
#' @param X
#' @param index
#'
#' @return
#'
#' @examples
q95 = function(X, index = 1:length(X)) {
  quantile(abs(X[index]), 0.95)
}
#' Title
#'
#' @param X
#' @param index
#'
#' @return
#'
#' @examples
q95hd = function(X, index = 1:length(X)){
  # Quantile estimate by Harrell & Davis 1982
  hd(abs(X[index]), 0.95)
}
#' Title
#'
#' @param X
#' @param index
#'
#' @return
#'
#' @examples
P1 = function(X, index = 1:length(X)) {
  sum(abs(X[index]) < eps) / length(index)
}
#' Title
#'
#' @param X
#' @param index
#'
#' @return
#'
#' @examples
W = function(X, index = 1:length(X)) {
  if (length(index) > 5000)
    index = sample(index, 5000)
  shapiro.test(X[index])$statistic
}
#' Title
#'
#' @param data
#' @param index
#' @param fscore
#' @param ...
#'
#' @return
#'
#' @examples
frank = function(data,index=1:nrow(data),fscore,...){
  S = apply(data[index,],2, fscore)
  order(S, decreasing = FALSE)
}
#' Title
#'
#' @param data
#' @param index
#' @param ...
#'
#' @return
#'
#' @examples
fRankMSIP = function(data, index=1:nrow(data), ...){
  nm  = ncol(data)
  N   = nrow(data)
  sip = matrix(NA, nrow = nm, ncol = nm)
  for (i in 1:nm) {
    vi = abs(data[index,i])
    for (j in 1:nm) {
      if (j==i) next
      vj = abs(data[index,j])
      sip[i, j] = mean( (vi - vj) < 0 )
    }
  }
  msip = rowMeans(sip, na.rm=TRUE)
  order(msip, decreasing = TRUE)
}
#' Title
#'
#' @param X
#' @param index
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fcor = function(X, index=1:length(X),...){
  cor(X[index,1],X[index,2])
}
#' Title
#'
#' @param X
#' @param index
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fcov = function(X, index=1:length(X),...){
  cov(X[index,1],X[index,2])
}
#' Title
#'
#' @param X
#' @param index
#' @param fscore
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fdif = function(X, index=1:nrow(X),fscore,...){
  fscore(X[,1],index,...) - fscore(X[,2],index,...)
}
#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
my5num = function(X) {
  c(
    hd(X, 0.05),
    hd(X, 0.25),
    hd(X, 0.5),
    hd(X, 0.75),
    hd(X, 0.95)
  )
}
