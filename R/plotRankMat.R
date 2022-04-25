#' Plots of Ranking Probability matrix Pr
#'
#' Interface to \code{corrplot::corrplot()} or mode/CI graph
#'
#' @param E
#' @param tab
#' @param score
#' @param type
#' @param method
#' @param nMC
#' @param cex.lab
#' @param show.main
#' @param offset
#' @param M
#' @param gPars
#'
#' @return
#' @export
#'
#' @examples
plotRankMat = function (
  E,
  tab = NULL,
  score='mue',
  type = 'levels',
  method = 'square',
  nMC = 1000,
  cex.lab = 1,
  show.main = TRUE,
  offset = 0.7,
  M = nrow(E),
  label  = 0,
  gPars = ErrViewLib::setgPars()
) {

  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(
    mfrow = c(1, 1),
    pty = pty,
    # mar = c(1,1,1,1),
    mgp = mgp,
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  if(is.null(tab)) {
    if( M == nrow(E) )
      tab = rankBS(E, score, nMC)
    else
      tab = rankBS2(E, score, nMC, M)
  }

  x = 1:ncol(E)
  main = paste0(toupper(score),' ranks (N = ',nrow(E),')')

  if(type == 'levels') {
    col2 <- colorRampPalette(
      c("#67001F", "#B2182B", "#D6604D", "#F4A582",
        "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
        "#4393C3", "#2166AC", "#053061"))
    colors = col2(200)[101:200]

    if(!show.main)
      mar = c(0,0,0,0) # else mar is provided by gPars

    corrplot::corrplot(
      tab$pRank[tab$mRank,],
      col.lim = c(0,1),
      col = colors,
      is.corr = FALSE,
      diag   = TRUE,
      method = method,
      order  = 'original',
      tl.col = 'black',
      tl.srt = 0,
      tl.offset = offset,
      tl.cex = cex.lab,
      cl.cex = cex.lab,
      win.asp = 1,
      mar = mar
    )
    if(show.main)
      mtext(main,
            line=0.70,
            cex = 1.5*cex.lab*cex,
            font = 2,
            adj = 0.2)
    if(label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj  = 0.975,
        # las  = 1,
        cex  = 1.5*cex.lab*cex,
        line = 0.70)

  } else {
    par(mar = c(3.5,7,3.5,1),
        pty='s',
        xaxt='n',
        yaxt='n')
    plot(x,x,type = 'n',xlab='',ylab='',main = '')
    grid()
    abline(v=1:ncol(E), lty=2, col= 'gray70')
    prank = tab$pRank[tab$mRank,]
    for(im in tab$mRank) {
      pt = cumsum(prank[im,])
      vmin = which(pt > 0.05)[1]
      vmax = which(pt > 0.95)[1]
      y = ncol(E)-im+1
      segments(vmin,y,vmax,y,lwd=2*lwd,col=cols[5])
      points(which.max(prank[im,]),y, pch = 18,
             cex=2, col = cols[1])
    }
    mtext(text = colnames(E)[tab$mRank], side = 2, adj= 1, line=0.3,
          at = rev(1:ncol(E)), las=2, cex = cex*cex.lab)
    mtext(text = 1:ncol(E), side = 3, line = 0.1*cex.lab,
          at = 1:ncol(E), las=1, cex = cex*cex.lab)
    box()
    if(show.main)
      mtext(main,
            line = 2.7,
            cex = 1.5*cex.lab*cex,
            font = 2,
            adj = 0.2)
    if(label > 0)
      mtext(
        text = paste0('(', letters[label], ')'),
        side = 3,
        adj  = 0.975,
        # las  = 1,
        cex  = 1.5*cex.lab*cex,
        line = 2.7)
  }

}