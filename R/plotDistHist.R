#' Scatter plot of an error set with marginal histogram
#'
#' @param x
#' @param y
#' @param uy
#' @param nclass
#' @param xlab
#' @param ylab
#' @param plotGauss
#' @param outLiers
#' @param p
#' @param labels
#' @param select
#' @param main
#' @param plotReg
#' @param plotConf
#' @param plotBA
#' @param plotBAci
#' @param xlim
#' @param ylim
#' @param scaleLegBA
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotDistHist = function(
  x,
  y,
  uy        = NULL,
  nclass    = NULL,  # Nb class for histogram
  xlab      = 'x',
  ylab      = 'y',
  plotGauss = FALSE,# Plot Gaussian fit of hist.
  outLiers  = FALSE, # Mark outliers
  p         = 0.9,   # Width of proba interval to detect outliers
  labels    = 1:length(x),
  select    = NULL,  # Indices of points to colorize
  main      = NULL,
  plotReg   = TRUE,  # Regression line
  plotConf  = FALSE, # Confidence limits on reg-line
  plotBA    = FALSE, # Bland-Altman LOAs
  plotBAci  = FALSE, # 95% CI on Bland-Altman LOAs
  xlim      = range(x),
  ylim      = range(y),
  scaleLegBA= 0.75,
  gPars
) {

  if (length(x)*length(y) == 0)
    return()

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mgp = mgp,
    pty = 'm',
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  colp = cols_tr2[5]
  if (!is.null(select)) {
    y1 = y[select]
    y2 = y[!select]
    colp = rep(cols_tr2[2], length(select))
    colp[!select] = cols_tr2[5]
  }

  # Subfigure with histogram
  par(mar = c(3, 3, 1.6, 0),
      fig = c(0, 0.35, 0, 1))
  h = hist(y, nclass = nclass, plot = FALSE)
  binWidth = h$breaks[2] - h$breaks[1]
  n = length(h$counts)
  x.l = rep(0, n)
  x.r = x.l - h$counts
  y.b = h$breaks[1:n]
  y.t = h$breaks[2:(n + 1)]
  plot(
    x.l,
    y.t,
    type = 'n',
    ylim = ylim,
    ylab = ylab,
    xlim = c(-1.1 * max(h$counts), 0),
    xaxt = 'n',
    xaxs = 'i',
    xlab = ''
  )
  grid()
  rect(
    xleft = x.l,
    ybottom = y.b,
    xright = x.r,
    ytop = y.t,
    border = NA,
    col = cols_tr2[5]
  )
  if (plotGauss) {
    ym = mean(y)
    ys = sd(y)
    xg  = seq(ym - 6 * ys, ym + 6 * ys, length.out = 1000)
    yg  = dnorm(xg, ym, ys)
    yg  = yg / max(yg) * max(h$counts)
    lines(-yg, xg, col = cols[6])
  }
  abline(h = 0, lty = 3)

  if (!is.null(select)) {
    y1 = y[select]
    h = hist(y1, breaks = h$breaks, plot = FALSE)
    n = length(h$counts)
    x.l = rep(0, n)
    x.r = x.l - h$counts
    y.b = h$breaks[1:n]
    y.t = h$breaks[2:(n + 1)]
    rect(
      xleft = x.l,
      ybottom = y.b,
      xright = x.r,
      ytop = y.t,
      density = -1,
      border = NA,
      col = cols_tr2[2]
    )

    y1 = y[!select]
    h = hist(y1, breaks = h$breaks, plot = FALSE)
    n = length(h$counts)
    x.l = rep(0, n)
    x.r = x.l - h$counts
    y.b = h$breaks[1:n]
    y.t = h$breaks[2:(n + 1)]
    rect(
      xleft = x.l,
      ybottom = y.b,
      xright = x.r,
      ytop = y.t,
      density = -1,
      border = NA,
      col = cols_tr2[5]
    )
  }
  box()

  par(
    mar = c(3, 0, 1.6, ifelse(plotBA,3,0.5)),
    fig = c(0.35, 1, 0, 1),
    new = TRUE
  )
  pch = 16

  # Transparent filled symbols
  if (!is.null(select)) {
    y1 = y[select]
    y2 = y[!select]
    colp = rep(cols_tr2[2], length(select))
    colp[!select] = cols_tr2[5]
  }
  plot(
    x,
    y,
    pch = pch,
    col = colp,
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    yaxt = 'n',
    cex = 0.75*cex,
    main = NULL
  )
  grid()
  if (!is.null(uy))
    segments(x, y - 2 * uy, x, y + 2 * uy, col = colp)
  nClass = length(unique(colp))
  legend('topright',
         title = main,
         bty = 'n',
         legend = '')
  abline(h = 0, lty = 3)

  if (outLiers) {
    # Mark and label quantile-based outliers
    plow = (1 - p) / 2
    pup  = p + plow
    lab  = y > quantile(y, p = pup) | y < quantile(y, p = plow)
    if (sum(lab) > 0) {
      points(
        x = x[lab],
        y = y[lab],
        pch = 16,
        col = cols[5]
      )
      text(
        x = x[lab],
        y = y[lab],
        labels = labels[lab],
        pos = 4
      )
    }
  }

  if (plotReg) {
    # Plot regression line
    reg = lm(y ~ x)
    indx = order(x)

    if(plotConf) {
      # Plot 95% confidence interval on reg-line
      p = predict(reg, interval = 'conf')
      matlines(x[indx],
               p[indx, ],
               col = cols[1],
               lwd = gPars$lwd,
               lty = c(1, 2, 2))
    } else {
      # Plot only regline
      p = predict(reg)
      matlines(x[indx],
               p[indx],
               col = cols[1],
               lwd = gPars$lwd,
               lty = 1)
    }

  }

  if(plotBA) {
    # Bland-Alman-type plot with LOAs
    bias = mean(y)
    abline(h = bias, col=cols[3])
    mtext(
      'Mean',
      side = 4,
      at = bias,
      col = cols[3],
      cex = scaleLegBA * cex,
      las = 1,
      line=0.25
    )

    if(plotBAci) {
      # 95% CI on mean
      ubias = sd(y)/sqrt(length(y))
      xlim = range(pretty(x))
      polygonXlimits = c(xlim, rev(range(xlim)))
      polygon(polygonXlimits,
              c(bias-1.96*ubias, bias-1.96*ubias,
                bias+1.96*ubias,bias+1.96*ubias),
              col = cols_tr2[3], border = NA)
    }

    # LOAs
    loas = quantile(y,probs = c(0.025,0.975))
    abline(h = loas, col=cols[c(2,4)])
    mtext(
      '02.5%',
      side = 4,
      at = loas[1],
      col = cols[2],
      cex =  scaleLegBA * cex,
      las = 1,
      line=0.25
    )
    mtext(
      '97.5%',
      side = 4,
      at = loas[2],
      col = cols[4],
      cex =  scaleLegBA * cex,
      las = 1,
      line=0.25
    )
    if(plotBAci) {
      # Bootstrap 95% CI on LOAs
      q = function(x,i) quantile(x[i],p=0.025)
      loas.boot = boot::boot(y, q, stype='i', R=1000)
      loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="basic")
      polygon(polygonXlimits,
              c(loas.ci$basic[4], loas.ci$basic[4],
                loas.ci$basic[5], loas.ci$basic[5]),
              col = cols_tr2[2], border = NA)

      q = function(x,i) quantile(x[i],p=0.975)
      loas.boot = boot::boot(y, q, stype='i', R=1000)
      loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="basic")
      polygon(polygonXlimits,
              c(loas.ci$basic[4], loas.ci$basic[4],
                loas.ci$basic[5], loas.ci$basic[5]),
              col = cols_tr2[4], border = NA)
    }
  }
  box()
}
