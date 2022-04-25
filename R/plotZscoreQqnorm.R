#' plotZscoreQqnorm
#'
#' @param R
#' @param sig
#' @param lim
#' @param title
#' @param gPars
#' @param plotCI
#' @param plotReg
#' @param markOutliers
#' @param labels
#'
#' @return
#' @export
#'
#' @examples
plotZscoreQqnorm = function(
  R,
  sig,
  lim = NULL,
  plotCI = FALSE,
  plotReg = FALSE,
  markOutliers = FALSE,
  labels = 1:length(R),
  title = '',
  gPars = ErrViewLib::setgPars()
) {

  if(plotCI) {
    # Monte Carlo estimation of 95% CI of qqnorm
    # for a given sample size N
    N   = length(R)
    nMC = 1000
    uly = matrix(0, ncol = N, nrow = nMC)
    for (i in 1:nMC) {
      t = sort(rnorm(N))
      qt = qqnorm(t, plot.it = FALSE)
      uly[i, ] = qt$y
    }
    q95 = t(apply(uly, 2, function(x)
      quantile(x, probs = c(0.025, 0.975))))
    rm(uly)
    ulx = sort(qt$x) # Set of quantiles for N points
  }

  # Z-scores
  Zs = R / sig

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2
  )

  if (is.null(lim))
    lim = max(3, max(abs(Zs)))

  q = qqnorm(
    Zs,
    main = '',
    pch = 16,
    col = cols[5],
    xlim = lim * c(-1, 1),
    ylim = lim * c(-1, 1)
  )
  grid()
  abline(a = 0, b = 1)

  if(plotCI)
    polygon(c(ulx, rev(ulx)),
            c(q95[, 1], rev(q95[, 2])),
            border = NA,
            col = cols_tr2[3])

  points(q$x,
         q$y,
         pch = 16,
         cex = 0.8,
         col = cols[5])

  if(plotReg)
    qqline(Zs, col = 2)

  if(markOutliers & plotCI) {
    out = (Zs - q95[, 1]) * (Zs - q95[, 2]) > 0
    text(x = ulx[out],
         y = Zs[out],
         labels = labels[out])
  }

  legend(
    'topleft',
    inset = 0.025,
    title = title,
    title.adj = 0,
    bty = 'n',
    legend = c('Data',
               '95% CI',
               'QQ line'),
    pch = c(16, -1, -1),
    col = c(cols[5], cols_tr2[3], 2),
    lty = c(NA, 1, 1),
    lwd = c(0, 10*lwd, lwd)
  )
  box()
}