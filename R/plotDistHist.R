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
#' @param colPoints
#' @param plotReg
#' @param plotConf
#' @param degree
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
  plotGauss = FALSE, # Plot Gaussian fit of hist.
  outLiers  = FALSE, # Mark outliers
  p         = 0.9,   # Width of proba interval to detect outliers
  labels    = 1:length(x),
  select    = NULL,  # Indices of points to colorize
  main      = NULL,
  plotReg   = TRUE,  # Regression line
  plotConf  = FALSE, # Confidence limits on reg-line
  degree    = 0,
  colPoints = NULL,
  plotBA    = FALSE, # Bland-Altman LOAs
  plotBAci  = FALSE, # 95% CI on Bland-Altman LOAs
  xlim      = range(x),
  ylim      = range(y),
  scaleLegBA = 0.75,
  scalePoints = 0.75,
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
    col = if(!is.null(colPoints)) colPoints else colp,
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    yaxt = 'n',
    cex = scalePoints*cex,
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
    lab  = y > ErrViewLib::hd(y, pup) | y < ErrViewLib::hd(y, plow)
    if (sum(lab) > 0) {
      points(
        x = x[lab],
        y = y[lab],
        pch = 16,
        col = cols[5],
        cex = scalePoints*cex
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
    # Build regression formula
    fo = y ~ 1
    if (degree > 0)
      fo = as.formula(
        paste0('y ~ 1 +',
               paste0(
                 'I(x^', 1:degree, ')',
                 collapse = '+'
               )))
    reg = lm(fo)
    indx = order(x)

    if(plotConf) {
      # Plot 95% confidence interval on reg-line
      p = predict(reg, interval = 'conf')
      matlines(x[indx],
               p[indx, ],
               col = cols[2],
               lwd = gPars$lwd,
               lty = c(1, 2, 2))
    } else {
      # Plot only regline
      p = predict(reg)
      matlines(x[indx],
               p[indx],
               col = cols[2],
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
                bias+1.96*ubias, bias+1.96*ubias),
              col = cols_tr2[3], border = NA)
    }

    # LOAs
    loas = c(
      ErrViewLib::hd(y, 0.025),
      ErrViewLib::hd(y, 0.975)
    )
    # loas = quantile(y,probs = c(0.025,0.975))
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
      q = function(x,i) ErrViewLib::hd(x[i], 0.025)
      # Adapt nBoot for boot.ci
      # (https://stat.ethz.ch/pipermail/r-help/2011-February/269006.html)
      nBoot = max(1000, length(y) + 1)
      loas.boot = boot::boot(y, q, stype='i', R=nBoot)
      loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="bca")
      polygon(polygonXlimits,
              c(loas.ci$bca[4], loas.ci$bca[4],
                loas.ci$bca[5], loas.ci$bca[5]),
              col = cols_tr2[2], border = NA)

      q = function(x,i) ErrViewLib::hd(x[i], 0.975)
      loas.boot = boot::boot(y, q, stype='i', R=nBoot)
      loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="bca")
      polygon(polygonXlimits,
              c(loas.ci$bca[4], loas.ci$bca[4],
                loas.ci$bca[5], loas.ci$bca[5]),
              col = cols_tr2[4], border = NA)
    }
  }
  box()
}
#' Scatter plot of an error set with marginal histogram (plotly version)
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
#' @param degree
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
plotlyDistHist = function(
  x,
  y,
  uy        = NULL,
  nclass    = NULL,  # Nb class for histogram
  xlab      = 'x',
  ylab      = 'y',
  plotGauss = FALSE, # Plot Gaussian fit of hist.
  outLiers  = FALSE, # Mark outliers
  p         = 0.9,   # Width of proba interval to detect outliers
  labels    = 1:length(x),
  select    = NULL,  # Indices of points to colorize
  main      = NULL,
  plotReg   = TRUE,  # Regression line
  plotConf  = FALSE, # Confidence limits on reg-line
  degree    = 0,
  plotBA    = FALSE, # Bland-Altman LOAs
  plotBAci  = FALSE, # 95% CI on Bland-Altman LOAs
  xlim      = range(x),
  ylim      = range(y),
  scaleLegBA = 0.75,
  scalePoints = 0.75,
  gPars
) {

  if (length(x)*length(y) == 0)
    return()

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  # Histogram(s)
  ## Use R's hist breaks for compatibility with Gaussian plot
  h = hist(y, breaks = nclass, plot = FALSE)

  hi <- plot_ly(showlegend = FALSE) %>%
    add_histogram(
      y = y,
      # nbinsy = nclass,
      ybins = list(
        start =min(h$breaks),
        end = max(h$breaks)+1,
        size = diff(h$breaks)[1]
      ),
      name = 'Density',
      hoverinfo = 'x',
      marker = list(color = cols[5]),
      alpha = 0.5
    )

  if (plotGauss) {
    ## Gaussian fit
    xg = seq(from = min(y),
             to   = max(y),
             length.out = 1000)
    yg = dnorm(xg, mean(y), sd(y))
    yg = yg / max(yg) * max(h$counts) # Rescale to hist's max
    hi <- hi %>%
      add_paths(
        x = yg, y = xg,
        name = 'Gaussian Fit',
        line = list(color = cols[3], width = lwd)
      )
  }

  if (!is.null(select)) {
    # Partial histograms: TBD !!!!
    # y1 = y[select]
    # y2 = y[!select]
  }

  # Scatterplot
  sc <- plot_ly(
    type = 'scatter',
    mode = 'markers',
    showlegend = FALSE) %>%
    add_trace(
      x = x, y = y,
      text = labels,
      name = 'Errors',
      hoverinfo = 'text',
      marker = list(size = 12*scalePoints, color = cols[5]),
      error_y = if(is.null(uy))
        list()
      else
        list(type="data",array=1.96*uy,color = cols[5],
             thickness = 0.75 * lwd, width=0)
    )


  if (outLiers) {
    # Mark quantile-based outliers
    plow = (1 - p) / 2
    pup  = p + plow
    out  = y > ErrViewLib::hd(y,pup) | y < ErrViewLib::hd(y,plow)
    if (sum(out) > 0) {
      sc <- sc %>%
        add_trace(
          x = x[out], y = y[out],
          text = labels[out],
          hoverinfo = 'text',
          marker = list(size = 12*scalePoints, color = cols[2]),
          error_y = if(is.null(uy))
            list()
          else
            list(type="data",array=1.96*uy[out],color = cols[2],
                 thickness = 0.75 * lwd, width=0)
        )
    }
  }

  if (plotReg) {
    # Build regression formula
    fo = y ~ 1
    if (degree > 0)
      fo = as.formula(
        paste0('y ~ 1 +',
               paste0(
                 'I(x^', 1:degree, ')',
                 collapse = '+'
               )))
    reg = lm(fo)
    indx = order(x)

    if(plotConf) {
      # Plot 95% confidence interval on reg-line
      p = predict(reg, interval = 'conf')
      sc <- sc %>%
        add_lines(
          x = x[indx], y = p[indx,1],
          line = list(color = cols[2], width = lwd)
        )%>%
        add_lines(
          x = x[indx], y = p[indx,2],
          line = list(color = cols[2], width = lwd, dash = 'dash')
        )%>%
        add_lines(
          x = x[indx], y = p[indx,3],
          line = list(color = cols[2], width = lwd, dash = 'dash')
        )
    } else {
      # Plot only regline
      p = predict(reg)
      sc <- sc %>%
        add_lines(
          x = x[indx], y = p[indx],
          line = list(color = cols[2], width = lwd)
        )
    }

  }

  if(plotBA) {
    # Bland-Alman-type plot with LOAs
    bias = mean(y)
    sc <- sc %>%
      add_lines(
        x = range(pretty(x)), y = c(bias,bias),
        line = list(color = cols[3], width = lwd)
      ) %>%
      plotly::layout(
        annotations = list(
          x = max(pretty(x)), y = bias,
          xanchor = 'left',
          yanchor = 'middle',
          text = 'Mean',
          font = list(family = 'Arial',
                      size = scaleLegBA *16,
                      color = cols[3]),
          showarrow = FALSE
        )
      )

    if(plotBAci) {
      # 95% CI on mean
      ubias = sd(y)/sqrt(length(y))
      xlim = range(pretty(x))
      sc <- sc %>%
        add_lines(
          x = xlim,
          y = c(bias - 1.96 * ubias, bias - 1.96 * ubias),
          line = list(color = 'transparent'),
          showlegend = FALSE
        ) %>%
        add_lines(
          x = xlim,
          y = c(bias + 1.96 * ubias, bias + 1.96 * ubias),
          fill = 'tonexty',
          fillcolor = cols_tr2[3],
          line = list(color = 'transparent'),
          showlegend = FALSE
        )
    }
    # LOAs
    loas = c(
      ErrViewLib::hd(y, 0.025),
      ErrViewLib::hd(y, 0.975)
    )
    # loas = quantile(y,probs = c(0.025,0.975))
    sc <- sc %>%
      add_lines(
        x = range(pretty(x)), y = c(loas[1],loas[1]),
        line = list(color = cols[2], width = lwd)
      ) %>%
      plotly::layout(
        annotations = list(
          x = max(pretty(x)), y = loas[1],
          xanchor = 'left',
          yanchor = 'middle',
          text = '02.5%',
          font = list(family = 'Arial',
                      size = scaleLegBA *16,
                      color = cols[2]),
          showarrow = FALSE
        )
      ) %>%
      add_lines(
        x = range(pretty(x)), y = c(loas[2],loas[2]),
        line = list(color = cols[4], width = lwd)
      ) %>%
      plotly::layout(
        annotations = list(
          x = max(pretty(x)), y = loas[2],
          xanchor = 'left',
          yanchor = 'middle',
          text = '97.5%',
          font = list(family = 'Arial',
                      size = scaleLegBA *16,
                      color = cols[4]),
          showarrow = FALSE
        )
      )

    if(plotBAci) {
      # Bootstrap 95% CI on LOAs
      q = function(x,i) ErrViewLib::hd(x[i], 0.025)
      # Adapt nBoot for boot.ci
      # (https://stat.ethz.ch/pipermail/r-help/2011-February/269006.html)
      nBoot = max(1000, length(y) + 1)
      loas.boot = boot::boot(y, q, stype='i', R=nBoot)
      loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="bca")
      sc <- sc %>%
        add_lines(
          x = xlim,
          y = c(loas.ci$bca[4], loas.ci$bca[4]),
          line = list(color = 'transparent'),
          showlegend = FALSE
        ) %>%
        add_lines(
          x = xlim,
          y = c(loas.ci$bca[5], loas.ci$bca[5]),
          fill = 'tonexty',
          fillcolor = cols_tr2[2],
          line = list(color = 'transparent'),
          showlegend = FALSE
        )

      q = function(x,i) ErrViewLib::hd(x[i], 0.975)
      loas.boot = boot::boot(y, q, stype='i', R=nBoot)
      loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="bca")
      sc <- sc %>%
        add_lines(
          x = xlim,
          y = c(loas.ci$bca[4], loas.ci$bca[4]),
          line = list(color = 'transparent'),
          showlegend = FALSE
        ) %>%
        add_lines(
          x = xlim,
          y = c(loas.ci$bca[5], loas.ci$bca[5]),
          fill = 'tonexty',
          fillcolor = cols_tr2[4],
          line = list(color = 'transparent'),
          showlegend = FALSE
        )
    }
  }

  # Display
  marg_plot <- subplot(
    hi,sc,
    nrows = 1,
    widths = c(.2,.7),
    margin = 0,
    shareY = TRUE
  ) %>%
    plotly::layout(
      xaxis = list(
        showline = TRUE, linewidth = 2,
        autorange = "reversed",
        color = 'gray70',
        showgrid = TRUE,
        mirror = TRUE
      ),
      yaxis = list(
        visible = TRUE,
        zeroline = TRUE, zerolinewidth = 2, zerolinecolor = '#AAA',
        showline = TRUE, linewidth = 2,
        side = "left",
        color = 'gray70',
        mirror = TRUE,
        showgrid = TRUE,
        title = ylab,
        range = range(pretty(ylim))
      ),
      xaxis2 = list(
        visible = TRUE,
        showline = TRUE, linewidth = 2,
        zeroline = TRUE,
        color = 'gray70',
        showgrid = TRUE,
        mirror = TRUE,
        title = xlab
      )
    ) %>%
    config(
      displaylogo = FALSE,
      modeBarButtonsToRemove =
        c("zoomIn2d", "zoomOut2d",
          "select2d","lasso2d",
          "autoScale2d",
          "toggleSpikelines",
          "hoverClosestCartesian",
          "hoverCompareCartesian")
    )
  marg_plot
}