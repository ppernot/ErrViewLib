#' plotPpnorm
#'
#' Deviation of ECDF from normality
#'
#' @param X
#' @param title
#' @param scale
#' @param plotCI
#' @param score
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotPpnorm <- function(
  X,
  title = '',
  scale  = FALSE,
  plotCI = FALSE,
  score  = TRUE,
  gPars
) {

  q = seq(-6, 6, length.out = 100)

  if(plotCI) {
    # Monte Carlo estimation of 95% CI
    N = length(q)
    M = length(X)
    nMC = 1000
    uly = matrix(0, ncol = N, nrow = nMC)
    for (i in 1:nMC)
      uly[i, ] = ecdf(rnorm(M))(q)
    q95 = t(apply(uly, 2, function(x)
      quantile(x, probs = c(0.025, 0.975))))
    rm(uly)
    ulx = pnorm(q)
  }

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = 's',
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  trapz = function (
    x,
    y
  ) {
    if(length(x) == 0)
      return(0)
    idx = 2:length(x)
    return(
      as.double((x[idx] - x[idx - 1]) %*%
                  (y[idx] + y[idx - 1]))/2
    )
  }

  if(scale)
    X = scale(X, center = TRUE, scale = TRUE)

  pt = c(0,pnorm(q),1)
  pe = c(0,ecdf(X)(q),1)
  plot(
    pt, pe,
    type = 'l',
    xlim = c(0,1), xlab = 'Expected cumulative distribution',
    ylim = c(0,1), ylab = 'Observed cumulative distribution',
    xaxs = 'i',
    yaxs = 'i',
    col  = cols[3],
    main = title
  )
  grid()
  abline(
    a=0,b=1,
    lty=2,
    col = cols[6]
  )
  polygon(
    c(pt,0),c(pe,0),
    col = cols_tr[6],
    border = NA
  )
  if(plotCI)
    matlines(
      ulx, q95,
      col = cols[6],
      lty = 3
    )
  box()
  if(score) {
    ma = trapz(pt,abs(pe-pt))
    text(
      0.6, 0.075,
      paste0(
        "Deviation area = ",
        signif(ma,2)
      ),
      cex=1
    )
  }
}
