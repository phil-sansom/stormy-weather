#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(ncdf4)

## Load source
source("src/datetime.R")
source("src/distance.R")
source("src/distances.R")
source("src/image.R")
source("src/invertlat.R")
source("src/lonflip.R")
source("src/tracecontour.R")


filename = "../cyclone_id/tests/output/cyclones_20000102T00.nc"
varid = "idclust"

threshold = 250


## Open file
nci = nc_open(filename)
lon = as.numeric(nci$dim$lon$vals)
lat = as.numeric(nci$dim$lat$vals)
nlon = length(lon)
nlat = length(lat)

times = ncvar_get(nci, "time")
time_units = ncatt_get(nci, "time", "units")$value
calendar = ncatt_get(nci, "time", "calendar")$value
facets = strsplit(time_units, " ")[[1]]
units = strsplit(time_units, " ")[[1]][1]
if (length(facets) == 3) {
  t0 = strptime(facets[3], format = "%Y-%m-%d", tz = "UTC")
} else if (length(facets) == 4) {
  t0 = strptime(paste(facets[3:4], collapse = " "),
                format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
} else if (length(facets) == 5) {
  t0 = strptime(paste(facets[3:5], collapse = " "),
                format = "%Y-%m-%d %H:%M:%OS %z", tz = "UTC")
}
t0 = format(t0, "%Y%m%dT%H%M%S")
dates = extract.dates(t0, times, calendar, units)
nt = length(times)

dd = distance(c(0,0), c(0,1))
dy = (lat[2] - lat[1])*dd

dny = floor(threshold/dy)
mx = which.min(abs(lon - 180))

## Initialize storage
output = array(NA, c(nlon,nlat,nt))

## Loop over times
for (t in 1:nt) {
  
  ## Load data
  input = ncvar_get(nci, varid, c(1,1,t), c(-1,-1,1))

  ## Extract object ids
  ids = sort(unique(as.numeric(input)))
  ids = ids[ids > 0 & !is.na(ids)]
  nids = length(ids)

  ## Loop over ids
  for (i in 1:nids) {
    
    points = which(input == ids[i], arr.ind = TRUE)
    output[cbind(points,t)] = ids[i]
    buffer = matrix(0, nlon, nlat)
    buffer[points] = 1
    
    indices = trace.contour(buffer)
    indices[,1] = (indices[,1] - 1) %% nlon + 1
    
    yjm1 = -180
    
    ## Loop over indices
    for (j in 1:nrow(indices)) {
      
      x = indices[j,1]
      y = indices[j,2]
      
      ## Create stencil
      if (y != yjm1) {
        
        y1 = min(y + dny,nlat)
        y0 = max(y - dny, 1)
        
        dx  = (lon[2] - lon[1])*min(cos(pi*lat[y1]/180),cos(pi*lat[y0]/180))*dd
        dnx = floor(threshold/dx)
        
        x1 = min(mx + dnx,nlon)
        x0 = max(mx - dnx,1)
        
        grid = as.matrix(expand.grid(x0:x1,y0:y1))
        dists = distances(c(180,lat[y]), cbind(lon[grid[,1]],lat[grid[,2]]))
        
        stencil = grid[dists <= threshold,]
        stencil[,1] = stencil[,1] - mx
        stencil[,2] = stencil[,2] - y
        
      }
      
      ## Use stencil to expand object
      stencilj = stencil
      stencilj[,1] = stencil[,1] + x
      stencilj[,2] = stencil[,2] + y
      stencilj = stencilj[1 <= stencilj[,2] & stencilj[,2] <= nlat,]
      stencilj[,1] = (stencilj[,1] -1) %% nlon + 1
      output[cbind(stencilj,t)] = ids[i]
      
      yjm1 = y
      
    } ## j
    
  } ## i
  
} ## t

image(lon, lat, output[,,1])
contour(lon, lat, input, add = TRUE, levels = 1:max(ids))
