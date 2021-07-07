library(geosphere)

distances = function(a, b, r = 6378137) {
  
  a[1] = a[1] %% 360
  if(a[1] > 180)
    a[1] = a[1] - 360
  b[,1] = b[,1] %% 360
  b[b[,1] > 180,1] = b[b[,1] > 180,1] - 360
  dist = distVincentySphere(a, b, r = r)/1000

  return(dist)
  
}
