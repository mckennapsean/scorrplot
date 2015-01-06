 

scorr.pretty.plot <- function(data, 
    colors = brewer.pal("Set1", n=9), 
    alpha=0.5,
    v1 = scorr.get.primary(),
    v2 = scorr.get.secondary(),
    bg ="white"
    ){
  
  nc <- ncol(data);
  X <- t(scale(t(data[,1:(nc-1)]), center=T, scale=T))/sqrt(ncol(data)-1)

  p1 <- X %*% v1
  p2 <- X %*% v2

  cols = paste(colors, as.hexmode(as.integer(alpha*255)), sep="")
  plot(c(-1.05,1.05), c(-1.05, 1.05),  col=bg, bty="n", xaxt="n", yaxt="n", 
       xlab=NA, ylab=NA, bg=bg)
  points(p1, p2, col=cols[as.numeric(data[, nc])], pch=19)
        
  draw.circle(x=0, y=0, radius=1, nv=100, border="darkgray", lwd=2)

} 
