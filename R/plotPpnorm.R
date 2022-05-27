#' plotPpnorm
#'
#' Deviation of ECDF from normality
#'
#' @param X -
#' @param title -
#' @param scale -
#' @param plotCI -
#' @param score -
#' @param dist -
#' @param shape -
#' @param label (integer) index of letter for subplot tag
#' @param gPars -
#'
#' @return a plot and an invisible list of metrics `misCal`, `misCalUp`
#'   and `calErr`
#' @export
#'
plotPpnorm <- function(
  X,
  title  = '',
  scale  = FALSE,
  plotCI = TRUE,
  score  = TRUE,
  dist   = c('norm','t'),
  shape  = 2,
  label  = 0,
  gPars  = ErrViewLib::setgPars()
) {

  dist    = match.arg(dist)

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
  pt = switch(
    dist,
    norm = c(0,normalp::pnormp(
      q * sqrt(shape^(2/shape)*gamma(3/shape)/gamma(1/shape)),
      p = shape),1),
    t    = c(0,pt(
      q * sqrt(shape / (shape-2)),
      df = shape),1)
  )

  if(score) {
    # Monte Carlo estimation of 95% CI
    N = length(q)
    M = length(X)
    nMC = 1000
    # uly = matrix(0, ncol = N, nrow = nMC)
    misCal = c()
    for (i in 1:nMC){
      Y = switch(
        dist,
        norm = normalp::rnormp(M, p = shape) /
          sqrt(shape^(2/shape)*gamma(3/shape)/gamma(1/shape)),
        t    = rt(M, df = shape) /
          sqrt(shape / (shape-2))
      )
      pe  = ecdf(Y)(q)
      misCal[i] = trapz(pt,abs(c(0,pe,1)-pt))
    }
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
    col = cols_tr[5],
    border = NA
  )
  misCalUp = NA
  if(plotCI) {
    lwr = upr = c()
    for(i in seq_along(pe)) {
      ci = DescTools::BinomCI(pe[i]*M, M, method = 'wilsoncc')
      lwr[i] = ci[,2]
      upr[i] = ci[,3]
    }
    misCalUp = trapz(pt,abs(upr-lwr)) / 2
    matlines(pt,cbind(lwr,upr),col = cols[3], lty = 2, lwd = lwd)
  }
  box()

  misCal = trapz(pt,abs(pe-pt))
  calErr = sqrt(sum((pt-pe)^2))
  if(score) {
    text(
      0.75, 0.1,
      paste0(
        'MisCal   = ',signif(misCal,2),'\n',
        ifelse(plotCI,paste0('MisCalUp = ',signif(misCalUp,2),'\n'),''),
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
