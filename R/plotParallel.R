#' Parallel plot of datasets with outliers detection features
#'
#' @param X
#' @param maxPoints
#' @param labels
#' @param lab.thresh
#' @param colors
#' @param rescale
#' @param scramble
#' @param outliers
#' @param ylab
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotParallel = function (X, maxPoints = nrow(X),
                         labels = NULL,
                         lab.thresh = 0,
                         colors = NULL,
                         rescale = TRUE,
                         scramble = FALSE,
                         outliers = FALSE,
                         outLabCex = 1,
                         xlim = c(1, ncol(X)),
                         ylim = NULL,
                         units = 'a.u.',
                         ylab = "Errors",
                         gPars) {
  # Driver for paraPlot

  ## Recast data to matrix
  if (class(X) == 'list') {
    n = names(X)
    X = as.matrix(as.data.frame(X))
    colnames(X) = n
  }

  ## Leave zero-variance colums out
  sdX = apply(X, 2, sd)
  X1 = X[, sdX != 0]

  ## Expose graphical params
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = c(8,3,0.2,3),
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2
  )

  ## Define color gradient
  if(is.null(colors))
    colors = genColors(rowMeans(X1))

  if(rescale) {
    X = apply(X1,2,scale)
    rownames(X) = rownames(X1)
    X1 = X
  }

  out = paraPlot(
    X1,
    col = colors,
    lwd = lwd,
    las = 2,
    var.label = labels,
    lab.thresh = lab.thresh,
    rescale = rescale,
    scramble = scramble,
    outliers = outliers,
    outLabCex = outLabCex,
    xlim = xlim,
    ylim = ylim,
    ylab = ifelse(
      rescale,
      'Centered-Scaled Errors',
      paste0('Errors [',units,']')
    ),
    cols_tr2 = cols_tr2
  )
  return(invisible(out))

}
#' Title
#'
#' @param x
#'
#' @return
#'
#' @examples
rescaleFun = function(x) {
  # Rescale to [-1,1]
  2*((x - min(x, na.rm = TRUE)) /
       (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))-0.5)
}
#' Title
#'
#' @param x
#' @param col
#' @param lty
#' @param pch
#' @param las
#' @param var.label
#' @param lab.thresh
#' @param rescale
#' @param scramble
#' @param outliers
#' @param ylab
#' @param cols_tr2
#' @param ...
#'
#' @return
#'
#' @examples
paraPlot = function (x,
                     col = 1,
                     lty = 1,
                     pch = 16,
                     las = las,
                     var.label = NULL,
                     lab.thresh = 0,
                     rescale = FALSE,
                     scramble = FALSE,
                     outliers = "no",
                     outLabCex = 1,
                     xlim = c(1, ncol(x)),
                     ylim = NULL,
                     ylab = "",
                     cols_tr2 = NULL,
                     ...) {
  # Parallel plot (adapted from MASS::parcoord)

  # Perturbation to horizontal positions
  rx = matrix(rep(1:ncol(x)),nrow=ncol(x),ncol=nrow(x))
  if(scramble)
    rx = rx + rnorm(length(rx),0,0.1)
# print(xlim)
  matplot(
    rx,
    t(x),
    type = "l",
    col  = col,
    lty  = lty,
    xlim = xlim,
    xlab = "",
    ylim = ylim,
    ylab = ylab,
    axes = FALSE,
    ...)

  axis(
    1,
    at = 1:ncol(x),
    labels = colnames(x),
    las = las)

  if(rescale) {
    axis(2,
      at = seq(-5, 5, by = 1),
      labels = seq(-5, 5, by = 1),
      pos = 1,
      las = las)
    for (i in 1:ncol(x))
      lines(c(i, i), c(-5, 5), col = "grey70")
    abline(h=-5:5, col = "grey90", lty=2)
  } else {
    ticks = pretty(as.matrix(x))
    axis(2,
      at = ticks,
      labels = ticks,
      pos = 1,
      las = las)
    for (i in 1:ncol(x))
      lines(c(i, i), range(ticks), col = "grey70")
    abline( h = ticks, col = "grey90", lty=2)
  }
  matpoints(
    rx,
    t(x),
    col = col,
    pch = pch,
    cex = 0.8)

  # Threshold-based Labels
  if(rescale) {
    if(!is.null(var.label)) {
      at = x[,ncol(x)]
      sel = abs(at) > lab.thresh
      if(sum(sel)>0) {
        at = at[sel]
        lab1 = var.label[sel]
        mtext(lab1,
              cex  = outLabCex,
              col  = 4,
              side = 4,
              las  = 2,
              line = -0.2,
              at   = at
        )
      }
      at = x[,1]
      sel = abs(at) > lab.thresh
      if(sum(sel)>0) {
        at = at[sel]
        lab2 = var.label[sel]
        mtext(lab2[sel],
              cex  = outLabCex,
              col  = 4,
              side = 2,
              las  = 2,
              line = -0.2,
              at   = at[sel]
        )
      }
    }
  }

  # Outliers zone
  sel = NULL
  if(outliers != "no") {
    if(outliers == "iqr") {
      qlim = t(apply(x,2,quantile,probs=c(0.25,0.75)))
      dq = qlim[,2]-qlim[,1]
      qlim[,1] = qlim[,1] - 1.5*dq
      qlim[,2] = qlim[,2] + 1.5*dq
    } else if (outliers == 'ci90') {
      qlim = t(apply(x,2,quantile,probs=c(0.05,0.95)))
    } else {
      qlim = t(apply(x,2,quantile,probs=c(0.025,0.975)))
    }
    polygon(
      c(1:ncol(x),rev(1:ncol(x))),
      c(qlim[,1],rev(qlim[,2])),
      border = NA, col=cols_tr2[4])
    if(!is.null(var.label)) {
      ql1  = matrix(qlim[,1],ncol=ncol(x),nrow=nrow(x),byrow = TRUE)
      ql2  = matrix(qlim[,2],ncol=ncol(x),nrow=nrow(x),byrow = TRUE)
      sel = rowSums(x < ql1 | x > ql2) == ncol(x)

      if(sum(sel)>0) {
        lab1 = var.label[sel]
        mtext(lab1,
              cex  = outLabCex,
              col  = 2,
              side = 4,
              las  = 2,
              line = -0.2,
              at   = x[sel,ncol(x)]
        )
        mtext(lab1,
              cex  = outLabCex,
              col  = 2,
              side = 2,
              adj  =  1,
              las  = 2,
              line = -0.2,
              at   = x[sel,1]
        )
      }
    }
  }
  return(invisible(sel))
}

