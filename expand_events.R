#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(geosphere)
library(ncdf4)

## Constants
radius = 6378137

## Parse arguments
args = commandArgs(TRUE)
threshold = as.numeric(args[1])
varid     = as.character(args[2])
infile    = as.character(args[3])
outfile   = as.character(args[4])

## Open file
nci = nc_open(infile)

## Extract dimensions
lon  = as.numeric(nci$dim$lon$vals)
lat  = as.numeric(nci$dim$lat$vals)
time = ncvar_get(nci, "time")
time.units = ncatt_get(nci, "time", "units"   )$value
calendar   = ncatt_get(nci, "time", "calendar")$value
nx = length(lon)
ny = length(lat)
nt = length(time)

## Assume regular grid
dx = lon[2] - lon[1]
dy = lat[2] - lat[1]

## Total x domain size
nxt = 360/dx

## Extract attributes
global.attributes = ncatt_get(nci, 0)

dk  = radius*2*pi/360/1000  ## 1 degree in km at equator
dky = dk*dy                 ## y grid distance in km
dny = floor(threshold/dky)  ## Number of grid boxes to search in y direction

## Create stencil library
stencils = list()
for (i in 1:ny) {
  
  y1 = min(i + dny,ny)
  y0 = max(i - dny, 1)
  
  dkx = dx*min(cos(pi*lat[y1]/180),cos(pi*lat[y0]/180))*dk
  dnx = floor(threshold/dkx)  ## Number of grid boxes to search in x direction
  
  if (nxt/2 < nx & nx < nxt) {
    dnx = min(dnx,nx)
  } else {
    dnx = min(dnx,floor(nxt/2))
  }
  x0 = -dnx
  x1 = +dnx
  lon1 = seq(x0*dx,x1*dx,dx)
  
  grid = as.matrix(expand.grid(x0:x1,y0:y1))
  dists = distVincentySphere(c(0,lat[i]), 
                             cbind(lon1[grid[,1]+dnx+1],lat[grid[,2]]))/1000
  
  stencils[[i]] = grid[dists <= threshold,]
  
} ## i

## Attributes
atts = global.attributes[! names(global.attributes) %in% "history"]
make_missing_value = nci$var[[varid]]$make_missing_value

## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude" )
time.dim = ncdim_def("time", time.units, time,
                     unlim = TRUE, calendar = calendar, longname = "Time")

## Define variables
object_nc = ncvar_def(varid, nci$var[[varid]]$units, 
                      list(lon.dim,lat.dim,time.dim),  
                      if (make_missing_value) nci$var[[varid]]$missval else NULL,
                      longname = nci$var[[varid]]$longname, 
                      prec = nci$var[[varid]]$prec,
                      compression = 5)

## Create netCDF file
nco = nc_create(outfile, list(object_nc))

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
                 threshold, " ", varid, " ", infile, " ", outfile)
if ("history" %in% names(global.attributes))
  history = paste(history, global.attributes$history, sep = "\n")
ncatt_put(nco, 0, "history", history)

## Loop over times
for (t in 1:nt) {
  
  ## Initialize storage
  output = matrix(0, nx, ny)
  
  ## Load data
  input = ncvar_get(nci, varid, c(1,1,t), c(-1,-1,1))
  
  ## Identify points for expansion
  indices = which(input > 0, arr.ind = TRUE)
  
  ## Loop over indices
  for (i in 1:nrow(indices)) {
    
    x = indices[i,1]
    y = indices[i,2]
    
    ## Use stencil to expand object
    stencil = stencils[[y]]
    stencil[,1] = stencil[,1] + x
    if (nx < nxt) {
      stencil = stencil[1 <= stencil[,1] & stencil[,1] <= nx,]
    } else {
      stencil[,1] = (stencil[,1] - 1) %% nx + 1
    }
    output[stencil] = 1
    
  } ## i
  
  ## Write data
  ncvar_put(nco, varid, output, c(1,1,t), c(-1,-1,1))
  
} ## t

## Close output file
nc_close(nco)

## Close input file
nc_close(nci)
