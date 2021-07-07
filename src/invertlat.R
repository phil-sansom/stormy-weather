invertlat = function(x) {
  
  y = array(NA, dim(x))
  nn = ncol(x)
  for (i in 1:nn) {
    y[,i] = x[,nn - i + 1]
  } ## i
  
  return(y)
  
}