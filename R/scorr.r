


scorr <- function (data, perms = 0, threshold = 0, coloring = 1:nrow(data),
    useDensity = T, showProfile = T, showPatch = F ){

  nc <- ncol(data)  

  v <- apply(data[,1:(nc-1)], 1, "var")
  if(sum(v==0) > 0){
    warning("Data contains rows with zero variance. Removing zero variance rows. Changed data frame will be returned.")
  }
  data = data[v > threshold, ]
  coloring = coloring[v > threshold  ]


  l <-  data[, nc]
  lnames <- levels(l)
  ml = length(lnames)
  Y = as.matrix(data[, 1:(nc-1)])
  Y[] = as.double(Y[])

  .Call("scorr", ncol(Y), nrow(Y), t(Y), 
                      as.integer(as.numeric(l)-1), as.character(lnames), 
                      as.integer(perms), as.numeric(coloring),
                      as.integer(useDensity), as.integer(showProfile),
                      as.integer(showPatch) )
   invisible(data)
}

scorr.close <-  function(){
 .Call("scorrClose")
}

scorr.highlight.index <- function(indices){
  .Call("scorrHighlightIndex", as.integer(indices),
as.integer(length(indices)) )
}

scorr.highlight.name <- function(names){
  .Call("scorrHighlightName", as.character(names),
as.integer(length(names)) )
}

scorr.get.size <- function(){
  .Call("scorrGetSize")
}

scorr.get.selected <- function(){
  .Call("scorrGetSelected")
}

scorr.get.density <- function(){
  n <- .Call("scorrGetSize")
  ind <- .Call("scorrGetCorIndex", n)
  d <- .Call("scorrGetDensity", as.integer(ind), as.integer(length(ind)) )
  name <- scorr.get.name(ind)
  val <- .Call("scorrGetCorValue", as.integer(ind), as.integer(length(ind)) )
  df <- data.frame(index = ind, density = d, name=name, cor=val)
}

scorr.get.corr <- function(p){
  n <- .Call("scorrGetSize")
  ind <- .Call("scorrGetCorIndex", n)
  val <- .Call("scorrGetCorValue", as.integer(ind), as.integer(length(ind)) )
  name <- scorr.get.name(ind)
  df <- data.frame(index = ind, cor=val, name=name)
  df <- subset(df, cor > p)
}

scorr.get.acorr <- function(p){
  n <- .Call("scorrGetSize")
  ind <- .Call("scorrGetCorIndex", -n)
  val <- .Call("scorrGetCorValue", as.integer(ind), as.integer(length(ind)) )
  name <- scorr.get.name(ind)
  df <- data.frame(index = ind, cor=val, name=name)
  df <- subset(df, cor < -p)
}

scorr.get.cor <- function(n){
  ind <- .Call("scorrGetCorIndex", as.integer(n) )
  val <- .Call("scorrGetCorValue", as.integer(ind), as.integer(length(ind)) )
  name <- scorr.get.name(ind)
  df <- data.frame(index = ind, cor=val, name=name)
}

scorr.get.acor <- function(n){
  ind <- .Call("scorrGetCorIndex", as.integer(-n) )
  val <- .Call("scorrGetCorValue", as.integer(ind), as.integer(length(ind)) )
  name <- scorr.get.name(ind);
  df <- data.frame(index = ind, cor=val, name=name)
}

scorr.get.name <- function(indices){
  .Call("scorrGetName", as.integer(indices), as.integer(length(indices)) );
}


scorr.set.primary <- function(v){
  .Call("scorrSetProjection", as.double(v), as.integer(length(v)), as.integer(0) )
}


scorr.set.secondary <- function(v){
  .Call("scorrSetProjection", as.double(v), as.integer(length(v)),
as.integer(1) )
}



scorr.get.primary <- function(v){
  .Call("scorrGetProjection", as.integer(0) )
}


scorr.get.secondary <- function(v){
  .Call("scorrGetProjection", as.integer(1) )
}

