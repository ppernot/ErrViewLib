#' Mean Signed Error (MSE)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
mse = function(X, index = 1:length(X),...) {
  mean(X[index])
}
#' Mode of signed errors by half-range estimator
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
hrmode = function(X, index = 1:length(X), ...) {
  genefilter::half.range.mode(X[index])
}
#' Mode of signed errors by half-sample estimator
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
hsmode = function(X, index = 1:length(X), ...) {
  modeest::hsm(X[index],bw=0.5)
}
#' Root-Mean Squared Deviation (RMSD)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
rmsd = function(X, index = 1:length(X),...) {
  sd(X[index])
}
#' Median Absolute Deviation (MAD-SD)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
mad_sd = function(X, index = 1:length(X),...) {
  mad(X[index])
}
#' Skewness (Skew)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
skew = function(X, index = 1:length(X),...) {
  moments::skewness(X[index])
}
#' Robust Skewness (SkewGM)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
skewgm = function(X, index = 1:length(X), ...) {
  X = X[index]
  m = hd(X, 0.5)
  (mean(X) - m) / mean(abs(X - m))
}
#' Kurtosis (Kurt)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
kurt = function(X, index = 1:length(X),...) {
  moments::kurtosis(X[index])
}
#' Robust Kurtosis (Crow & Siddiqui)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
kurtcs = function(X, index = 1:length(X), ...) {
  # Formula KR4 from Kim2004
  x = X[index]
  q975 = hd(x, 0.975)
  q025 = hd(x, 0.025)
  q75  = hd(x, 0.75)
  q25  = hd(x, 0.25)

  (q975-q025)/(q75-q25) - 2.91
}
#' Mean Unsigned Error (MUE)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
mue = function(X, index = 1:length(X),...) {
  mean(abs(X[index]))
}
#' 95th Quantile of absolute errors (Q95)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
q95 = function(X, index = 1:length(X),...) {
  quantile(abs(X[index]), 0.95, na.rm = TRUE)
}
#' 95th Quantile of absolute errors using the HD algorithm (Q95HD)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
q95hd = function(X, index = 1:length(X),...){
  # Quantile estimate by Harrell & Davis 1982
  hd(abs(X[index]), 0.95)
}
#' Gini coefficient of absolute errors
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
gini = function(X, index = 1:length(X),...){
  ineq::Gini(abs(X[index]))
}
#' Gini coefficient of mode-centered absolute errors
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
gimc = function(X, index = 1:length(X), ...) {
  X = X[index]
  ineq::Gini(abs(X - hrmode(X)))
}

#' Pietra's coefficient of absolute errors
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
pietra = function(X, index = 1:length(X),...){
  ineq::RS(abs(X[index]))
}
#' Lorenz curve asymmetry of absolute errors
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
lasym = function(X, index = 1:length(X),...){
  ineq::Lasym(abs(X[index]))
}
#' Zanardi index (Clementi2019)
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
Zanardi = function(X, index = 1:length(X), ...) {
  # Estimate Zanardi index (Clementi2019)

  trapz = function(x,y) {
    # Trapezoidal integration
    idx = 2:length(x)
    return (as.double( (x[idx] - x[idx-1]) %*%
                         (y[idx] + y[idx-1])) / 2)
  }
  lagint = function(x,y,xout) {
    # Lagrange quadratic interpolation
    x1 = x[1]; y1 = y[1]
    x2 = x[2]; y2 = y[2]
    x3 = x[3]; y3 = y[3]
    yout = (xout-x2)/(x1-x2) * (xout-x3)/(x1-x3) * y1 +
      (xout-x1)/(x2-x1) * (xout-x3)/(x2-x3) * y2 +
      (xout-x1)/(x3-x1) * (xout-x2)/(x3-x2) * y3
    return(yout)
  }

  X = abs(X[index])

  # Lorenz curve
  L = ineq::Lc(X, plot = FALSE)

  # Gini coefficient
  G = ineq::Gini(X)

  # Gini transform
  p = L$p
  q = L$L

  x = p + q
  y = p - q

  # Discriminant point
  d = which.min(abs(x-1.0))
  xd = 1.0
  yd = lagint(x[(d-1):(d+1)],y[(d-1):(d+1)],xd)
  pd = 0.5 * (xd + yd)
  qd = 0.5 * (xd - yd)

  # Zanardi index
  n   = length(x)
  Kd  = pd * qd / 2

  i0  = d # Position of integration limit vs. d
  if(x[d]-xd > 0) i0 = d-1
  A0p = trapz(x[1:i0], y[1:i0]) - yd / 2
  A0r = trapz(x[(i0+1):n], y[(i0+1):n]) - yd / 2
  Gp  = A0p / Kd
  Gr  = A0r / Kd
  Zd  = 2 * Kd * (Gr - Gp) / G

  return(Zd)
}
#' Probability for the absolute errors to be below a threshold
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
P1 = function(X, index = 1:length(X), eps) {
  mean(abs(X[index]) < eps)
}
#' Normality index of a sample by the Shapiro test
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#'
#' @return
#' @export
#'
W = function(X, index = 1:length(X),...) {
  if (length(index) > 5000)
    index = sample(index, 5000)
  shapiro.test(X[index])$statistic
}
#' Order of a series of statistics
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param data -
#' @param index -
#' @param fscore -
#' @param ... -
#'
#' @return
#' @export
#'
forder = function(data,index=1:nrow(data),fscore,...){
  S = apply(data[index,],2, fscore)
  # Rq: might use rank instead of order, but the pRank matrix
  # would need to be transposed in rankBS...
  order(S, decreasing = FALSE)
}
#' Order of a series of MSIP statistics
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param data -
#' @param index -
#' @param ...
#'
#' @return
#' @export
#'
fOrderMSIP = function(data, index=1:nrow(data), ...){
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
#' Correlation
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#' @param ...
#'
#' @return
#' @export
#'
fcor = function(X, index=1:length(X),...){
  cor(X[index,1],X[index,2])
}
#' Covariance
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#' @param ... -
#'
#' @return
#' @export
#'
fcov = function(X, index=1:length(X),...){
  cov(X[index,1],X[index,2])
}
#' Difference of statistics
#'
#' Auxiliary function for bootstrap by 'boot::boot()'
#'
#' @param X -
#' @param index -
#' @param fscore -
#' @param ... -
#'
#' @return
#' @export
#'
fdif = function(X, index=1:nrow(X),fscore,...){
  fscore(X[,1],index,...) - fscore(X[,2],index,...)
}
#' 5-numbers summary of sample
#'
#' @param X -
#'
#' @return
#' @export
#'
my5num = function(X) {
  c(
    hd(X, 0.05),
    hd(X, 0.25),
    hd(X, 0.5),
    hd(X, 0.75),
    hd(X, 0.95)
  )
}
#' Mean of the Folded Normal distribution
#'
#' @param x -
#'
#' @return
#' @export
#'
muF = function(x) {
  mu=x[1]; sig=x[2]
  sig*sqrt(2/pi)*exp(-mu^2/(2*sig^2)) -mu*erf(-mu/(sqrt(2)*sig))
}
#' CDF of the Folded Normal distribution
#'
#' @param x -
#'
#' @return
#' @export
#'
cdfF = function(x) {
  u = x[1]; mu = x[2]; sig = x[3]
  (erf((u+mu)/(sqrt(2)*sig))+erf((u-mu)/(sqrt(2)*sig)))/2
}
#' Q95 statistics of the Folded Normal distribution
#'
#' @param x -
#'
#' @return
#' @export
#'
q95F = function(x) {
  mu=x[1]; sig=x[2]
  fz = function(x,mu,sig,prob) {
    cdfF(c(x,mu,sig)) - prob
  }
  mueF = muF(c(mu,sig))
  uniroot(f = fz, interval=c(mueF,mueF+6*sig),check.conv = TRUE,
          mu = mu, sig = sig, prob = 0.95)$root
}
#' Agresti-Coull estimation of uncertainty on CDF
#'
#' @param X -
#' @param n -
#'
#' @return
#' @export
#'
agrestiCoull = function(X,n) {
  p=(X+1/2)/(n+1)
  return(sqrt(p*(1-p)/(n+1)))
}
#' Erf function
#'
#' @param x -
#'
#' @return
#' @export
#'
erf = function(x) {
  2 * pnorm(x * sqrt(2)) - 1
}
