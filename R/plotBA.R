#' Bland-Altman-type plot of error samples
#'
#' @param data1
#' @param data2
#' @param ylim
#' @param title
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotBA = function(
  data1,
  data2,
  ylim = NULL,
  title = '',
  gPars = ErrViewLib::setgPars()
) {

  # Bootstapped version of the Bland-Altman graph
  delta = data2 - data1
  meand = (data1 + data2) / 2

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2
  )

  plot( meand, delta,
        pch=16, cex=0.5, col = cols[5],
        ylim =
          if(is.null(ylim))
            range(delta)
        else
          ylim,
        pty = 's',
        xlab = 'Means',
        ylab = 'Differences',
        main = title
  )
  grid()
  abline(h = 0)

  # Bias
  bias = mean(delta)
  abline(h = bias, col=cols[3])

  ubias = sd(delta)/sqrt(length(delta))
  xlim = range(pretty(meand))
  polygonXlimits = c(xlim, rev(range(xlim)))
  polygon(polygonXlimits,
          c(bias-1.96*ubias, bias-1.96*ubias,
            bias+1.96*ubias,bias+1.96*ubias),
          col = cols_tr2[3], border = NA)
  mtext(
    'Bias',
    side = 4,
    at = bias,
    col = cols[3],
    cex = 0.6 * cex,
    las = 1,
    line=0.25
  )

  # LOAs
  loas = quantile(delta,probs = c(0.025,0.975))
  abline(h = loas, col=cols[c(2,4)])
  mtext(
    '2.5%',
    side = 4,
    at = loas[1],
    col = cols[2],
    cex = 0.6 * cex,
    las = 1,
    line=0.25
  )
  mtext(
    '97.5%',
    side = 4,
    at = loas[2],
    col = cols[4],
    cex = 0.6 * cex,
    las = 1,
    line=0.25
  )

  q = function(x,i) ErrViewLib::hd(x[i],p=0.025)
  nBoot = max(1000, length(delta) + 1) # Needed for boot.ci
  loas.boot = boot::boot(delta, q, stype = 'i', R = nBoot)
  loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="bca")
  polygon(polygonXlimits,
          c(loas.ci$basic[4], loas.ci$basic[4],
            loas.ci$basic[5], loas.ci$basic[5]),
          col = cols_tr2[2], border = NA)

  q = function(x,i) ErrViewLib::hd(x[i],p=0.975)
  loas.boot = boot::boot(delta, q, stype = 'i', R = nBoot)
  loas.ci   = boot::boot.ci(loas.boot, conf = 0.95, type = "bca")
  polygon(polygonXlimits,
          c(loas.ci$basic[4], loas.ci$basic[4],
            loas.ci$basic[5], loas.ci$basic[5]),
          col = cols_tr2[4], border = NA)

}
