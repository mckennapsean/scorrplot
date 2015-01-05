
scorrplot <- function(X, primary, secondary, spearman=FALSE, c.lwd=2, c.col="black",
    g.lwd=1, g.locations = seq(-1,1, length.out=20), g.col="lightgray",
    xlab="Primary", ylab="", bty="n", pch=19, col="#00000050", ...){
  library(plotrix)

  if(spearman){
    X <- apply(X, 2, "order")
  }
  X <- scale( X ) / sqrt(nrow(X)-1)

  primary = primary - mean(primary)
  primary = primary / sqrt( sum(primary^2) )
  secondary = secondary - mean(secondary)
  secondary = secondary -  as.double((t(secondary)%*%primary )) * primary
  secondary = secondary/ sqrt(sum(secondary^2))


  U = cbind(primary, secondary)

  Xp = t(U) %*% X

  plot(NULL, xlim = c(-1, 1), ylim=c(-1, 1), xaxp=c(min(g.locations),
        max(g.locations), length(g.locations)),yaxp=c(min(g.locations),
        max(g.locations), length(g.locations)), asp=1, bty=bty, xlab=xlab,
        ylab=ylab, ...  )
  

  y = sqrt(1-g.locations^2)
  segments(x0 = g.locations, y0 = y, x1=g.locations,
      y1= -y, col=g.col, lwd=g.lwd)
  
  draw.circle(0, 0, 1, border=c.col, lwd=c.lwd)

  points( t(Xp), pch=pch, col=col, ...)


}

