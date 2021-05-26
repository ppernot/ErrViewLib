#' plotPpnorm
#'
#' Deviation of ECDF from normality
#'
#' @param X
#' @param title
#' @param scale
#' @param plotCI
#' @param score
#' @param label (integer) index of letter for subplot tag
#' @param gPars
#'
#' @return a plot and an invisible list of metrics `misCal`, `misCalUp`
#'   and `calErr`
#' @export
#'
#' @examples
plotPpnorm <- function(
  X,
  title = '',
  scale  = FALSE,
  plotCI = FALSE,
  score  = TRUE,
  label  = 0,
  gPars
) {

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

  q = seq(-6, 6, length.out = 100)
  pt = c(0,pnorm(q),1)

  if(plotCI) {
    # Monte Carlo estimation of 95% CI
    N = length(q)
    M = length(X)
    nMC = 1000
    uly = matrix(0, ncol = N, nrow = nMC)
    misCal = c()
    for (i in 1:nMC){
      pe  = ecdf(rnorm(M))(q)
      uly[i, ]  = pe
      misCal[i] = trapz(pt,abs(c(0,pe,1)-pt))
    }
    q95 = t(apply(uly, 2, function(x) ErrViewLib::vhd(x)))
    rm(uly)
    ulx = pnorm(q)
    misCalUp = ErrViewLib::hd(misCal, q = 0.975)
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

  if(scale)
    X = scale(X, center = TRUE, scale = TRUE)

  pe = c(0,ecdf(X)(q),1)
  plot(
    pt, pe,
    type = 'l',
    xlim = c(0,1), xlab = 'Expected cumulative distribution',
    ylim = c(0,1), ylab = 'Observed cumulative distribution',
    xaxs = 'i',
    yaxs = 'i',
    col  = cols[3],
    main = title,
    lwd  = lwd
  )
  grid()
  abline(
    a=0,b=1,
    lty=2,
    col = cols[6],
    lwd = lwd
  )
  polygon(
    c(pt,0),c(pe,0),
    col = cols_tr2[6],
    border = NA
  )
  if(plotCI)
    matlines(
      ulx, q95,
      col = cols[6],
      lwd = lwd,
      lty = 3
    )
  box()
  if(score) {
    misCal = trapz(pt,abs(pe-pt))
    calErr = sum((pt-pe)^2)
    text(
      0.75, 0.1,
      paste0(
        'MisCal   = ',signif(misCal,2),'\n',
        'MisCalUp = ',signif(misCalUp,2),'\n',
        'CalErr   = ',signif(calErr,2)
      ),
      cex=1
    )
  }

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)

  invisible(
    list(
      misCal = misCal,
      misCalUp = misCalUp,
      calErr = calErr
    )
  )
}
