# produce static s-CorrPlot for a dataset
# input data expects rows as variables & columns as observations
# by default, expects the two variables for the projection
# can also accept a column of factors for coloring the plot
scorr.plot <- function(
    data,
    primary,
    secondary,
    alpha = 0.1,
    ...){
  
  # load dependencies for plotting
  library(plotrix)
  library(RColorBrewer)
  
  # plotting parameters
  c.lwd = 2
  c.col = "black"
  g.lwd = 1
  g.lwd2 = 1.5
  g.locations = seq(-1, 1, length.out = 21)
  g.col = "lightgray"
  g.col2 = "darkgray"
  xlab = "correlation to variable of interest"
  ylab = ""
  bty = "n"
  pch = 19
  col = "#000000"
  col = paste(col, as.hexmode(as.integer(alpha * 255)), sep = "")
  
  # check input data for labels (factor type)
  # if present, these are used to color the variables
  # colors are automatically pulled from ColorBrewer library
  last.column.type = sapply(data, class)[ncol(data)]
  names(last.column.type) = NULL
  labels.exist = last.column.type == "factor"
  if(labels.exist){
    bg = "white"
    colors = brewer.pal("Set1", n = 9)
    cols = paste(colors, as.hexmode(as.integer(alpha * 255)), sep = "")
  }
  
  # determine our input data matrix, X
  # will vary depending on the column of factors
  if(labels.exist){
    X = t(data[, 1 : ncol(data) - 1])
  }else{
    X = t(data)
  }
  
  # standardize data matrix
  X = t(scale(X)) / sqrt(nrow(X) - 1)
  
  # for variables of interest, transpose when given row data
  # remove factor label from the variables if present
  if(!is.null(ncol(primary))){
    if(ncol(X) != ncol(primary)){
      primary = primary[, 1 : ncol(primary) - 1]
    }
    primary = t(primary)
  }
  if(!is.null(ncol(secondary))){
    if(ncol(X) != ncol(secondary)){
      secondary = secondary[, 1 : ncol(secondary) - 1]
    }
    secondary = t(secondary)
  }
  
  # standardize our primary and secondary variables of interest
  # secondary becomes the s-CorrPlot's y-axis projection vector
  primary = primary - mean(primary)
  primary = primary / sqrt( sum(primary^2) )
  secondary = secondary - mean(secondary)
  secondary = secondary -  as.double((t(secondary)%*%primary )) * primary
  secondary = secondary/ sqrt(sum(secondary^2))
  
  # projection plane (both vectors)
  # secondary is negative so that the secondary variable of
  # interest is on the bottom half of the plot
  U = cbind(primary, -secondary)
  
  # compute dot product of our projection plane & all the data
  Xp = t(U) %*% t(X)
  
  # prepare the plot
  plot(NULL, xlim = c(-1, 1), ylim=c(-1, 1), asp = 1, bty = bty, xlab = xlab,
        ylab = ylab, axes = FALSE, ...  )
  
  # axis for primary variable of interest
  axis(1, at = g.locations)
  
  # add gridlines
  y = sqrt(1 - g.locations ^ 2)
  segments(x0 = g.locations, y0 = y, x1 = g.locations,
           y1 = -y, col = g.col, lwd = g.lwd)
  segments(x0 = 0, y0 = 1, x1 = 0, y1 = -1, col = g.col2, lwd = g.lwd2)
  segments(x0 = -1, y0 = 0, x1 = 1, y1 = 0, col = g.col2, lwd = g.lwd2)
  
  # draw the s-CorrPlot projection plane circle
  draw.circle(0, 0, 1, border = c.col, lwd = c.lwd)
  
  # draw variables as points, color them if we were given factor labels
  if(labels.exist){
    points(t(Xp), col = cols[as.numeric(data[, ncol(data)])], pch = pch, ...)
  }else{
    points(t(Xp), pch = pch, col = col, ...)
  }
}
