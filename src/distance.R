library(geosphere)

distance = function(x, y, r = 6378137) {
  
  if (length(x) > 1) {
    xx = x %% 360
    xx[xx > 180] = xx[xx > 180] - 360
    dist = distVincentySphere(cbind(xx,y), r = r)/1000
  } else {
    dist = 0
  }
  return(dist)
  
}
