


gyroscope <- function (data, perms = 0, threshold = 0, coloring = 1:nrow(data),
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

  .Call("gyroscope", ncol(Y), nrow(Y), t(Y), 
                      as.integer(as.numeric(l)-1), as.character(lnames), 
                      as.integer(perms), as.numeric(coloring),
                      as.integer(useDensity), as.integer(showProfile),
                      as.integer(showPatch) )
   invisible(data)
}

gyroscope.close <-  function(){
 .Call("gyroscopeClose")
}

gyroscope.highlight.index <- function(indices){
  .Call("gyroscopeHighlightIndex", as.integer(indices),
as.integer(length(indices)) )
}

gyroscope.highlight.name <- function(names){
  .Call("gyroscopeHighlightName", as.character(names),
as.integer(length(names)) )
}

gyroscope.get.size <- function(){
  .Call("gyroscopeGetSize")
}

gyroscope.get.selected <- function(){
  .Call("gyroscopeGetSelected")
}

gyroscope.get.density <- function(){
  n <- .Call("gyroscopeGetSize")
  ind <- .Call("gyroscopeGetCorIndex", n)
  d <- .Call("gyroscopeGetDensity", as.integer(ind), as.integer(length(ind)) )
  name <- gyroscope.get.name(ind)
  val <- .Call("gyroscopeGetCorValue", as.integer(ind), as.integer(length(ind)) )
  df <- data.frame(index = ind, density = d, name=name, cor=val)
}

gyroscope.get.corr <- function(p){
  n <- .Call("gyroscopeGetSize")
  ind <- .Call("gyroscopeGetCorIndex", n)
  val <- .Call("gyroscopeGetCorValue", as.integer(ind), as.integer(length(ind)) )
  name <- gyroscope.get.name(ind)
  df <- data.frame(index = ind, cor=val, name=name)
  df <- subset(df, cor > p)
}

gyroscope.get.acorr <- function(p){
  n <- .Call("gyroscopeGetSize")
  ind <- .Call("gyroscopeGetCorIndex", -n)
  val <- .Call("gyroscopeGetCorValue", as.integer(ind), as.integer(length(ind)) )
  name <- gyroscope.get.name(ind)
  df <- data.frame(index = ind, cor=val, name=name)
  df <- subset(df, cor < -p)
}

gyroscope.get.cor <- function(n){
  ind <- .Call("gyroscopeGetCorIndex", as.integer(n) )
  val <- .Call("gyroscopeGetCorValue", as.integer(ind), as.integer(length(ind)) )
  name <- gyroscope.get.name(ind)
  df <- data.frame(index = ind, cor=val, name=name)
}

gyroscope.get.acor <- function(n){
  ind <- .Call("gyroscopeGetCorIndex", as.integer(-n) )
  val <- .Call("gyroscopeGetCorValue", as.integer(ind), as.integer(length(ind)) )
  name <- gyroscope.get.name(ind);
  df <- data.frame(index = ind, cor=val, name=name)
}

gyroscope.get.name <- function(indices){
  .Call("gyroscopeGetName", as.integer(indices), as.integer(length(indices)) );
}


gyroscope.set.primary <- function(v){
  .Call("gyroscopeSetProjection", as.double(v), as.integer(length(v)), as.integer(0) )
}


gyroscope.set.secondary <- function(v){
  .Call("gyroscopeSetProjection", as.double(v), as.integer(length(v)),
as.integer(1) )
}



gyroscope.get.primary <- function(v){
  .Call("gyroscopeGetProjection", as.integer(0) )
}


gyroscope.get.secondary <- function(v){
  .Call("gyroscopeGetProjection", as.integer(1) )
}

