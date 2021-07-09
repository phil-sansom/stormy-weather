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

## Parse arguments
args = commandArgs(TRUE)
filename  = as.character(args[1])
varid     = as.character(args[2])
threshold = as.numeric(args[3])

## Open file
nci = nc_open(filename)

## Extract dimensions
lon = as.numeric(nci$dim$lon$vals)
lat = as.numeric(nci$dim$lat$vals)
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
nlon = length(lon)
nlat = length(lat)
nt = length(times)

## Extract attributes
global.attributes = ncatt_get(nci, 0)

dd = distance(c(0,0), c(0,1))
dy = (lat[2] - lat[1])*dd

dny = floor(threshold/dy)
mx = which.min(abs(lon - 180))

## Initialize storage
output = array(0, c(nlon,nlat,nt))

## Loop over times
for (t in 1:nt) {
  
  ## Echo time
  print(t)
  
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


#############################
## Write expanded cyclones ##
#############################

## Attributes
atts = global.attributes[! names(global.attributes) %in% "history"]
make_missing_value = nci$var[[varid]]$make_missing_value

## Filename
outfile = tools::file_path_sans_ext(filename)
outfile = paste0(outfile, "_expanded.nc")

## Define dimensions
lon_nc  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat_nc  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude" )
time_nc = ncdim_def("time", time_units, times,
                    unlim = TRUE, calendar = calendar, longname = "Time")

## Define variables
object_nc = ncvar_def(varid, nci$var[[varid]]$units, 
                      list(lon_nc,lat_nc,time_nc),  
                      if (make_missing_value) nci$var[[varid]]$missval else NULL,
                      longname = nci$var[[varid]]$longname, 
                      prec = nci$var[[varid]]$prec,
                      compression = 5)

## Create netCDF file
nco = nc_create(outfile, list(object_nc))

## Write data
ncvar_put(nco, varid, output)

## Write standard names
ncatt_put(nco, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(nco, "latitude", "standard_name", "latitude", prec = "text")
ncatt_put(nco, "time", "standard_name", "time", prec = "text")

## Write global attributes
for (att in names(atts))
  ncatt_put(nco, 0, att, atts[[att]])
ncatt_put(nco, 0, "threshold", threshold)

## Write history
history = paste0(format(Sys.time(), "%FT%XZ%z"), ": ", "./expand_events.R ",
                 filename, " ", varid, " ", threshold)
if ("history" %in% names(global.attributes))
  history = paste(global.attributes$history, history, sep = "\n")
ncatt_put(nco, 0, "history", history)

## Close output file
nc_close(nco)

## Close input file
nc_close(nci)

