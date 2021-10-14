lonflip = function(x, lon) {

  if (!is.matrix(x))
    dim(x) = c(length(x),1)
  if (lon[1] < 0) {
    mask0to180 =    0 <= lon & lon < 180
    mask180to0 = -180 <= lon & lon <   0
    y = rbind(x[mask0to180,],x[mask180to0,])
    z = c(lon[mask0to180],lon[mask180to0] + 360)
  } else {
    mask0to180 =   0 <= lon & lon < 180
    mask180to0 = 180 <= lon & lon < 360
    y = rbind(x[mask180to0,],x[mask0to180,])
    z = c(lon[mask180to0] - 360,lon[mask0to180])
  }
  
  return(list(x = y, lon = z))
  
}
