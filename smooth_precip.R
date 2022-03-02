#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(geosphere)
library(ncdf4)

## Constants
radius = 6378137

## Parse arguments
args = commandArgs(TRUE)
threshold = as.numeric  (args[1L])
dinfile   = as.character(args[2L])
pinfile   = as.character(args[3L])
doutfile  = as.character(args[4L])
poutfile  = as.character(args[5L])

## Open files
ncid = nc_open(dinfile)
ncip = nc_open(pinfile)
vard = names(ncid$var)[1L]
varp = names(ncip$var)[1L]

## Extract dimensions
lon  = as.numeric(ncid$dim$lon$vals)
lat  = as.numeric(ncid$dim$lat$vals)
time = ncvar_get(ncid, "time")
time.units = ncatt_get(ncid, "time", "units"   )$value
calendar   = ncatt_get(ncid, "time", "calendar")$value
nx = length(lon)
ny = length(lat)
nt = length(time)

## Check dimensions
lonp  = as.numeric(ncip$dim$lon$vals)
latp  = as.numeric(ncip$dim$lat$vals)
timep = ncvar_get(ncip, "time")
if (!all(lonp == lon))
  stop (paste("Longitudes don't match beteween", dinfile, "and", pinfile))
if (!all(latp == lat))
  stop (paste("Latitudes don't match beteween", dinfile, "and", pinfile))
if (!all(timep == time))
  stop (paste("Times don't match beteween", dinfile, "and", pinfile))


##############
## Stencils ##
##############

## Assume regular grid
dx = lon[2L] - lon[1L]
dy = lat[2L] - lat[1L]

## Total x domain size
nxt = 360/dx

dk  = radius*2*pi/360/1000  ## 1 degree in km at equator
dky = dk*dy                 ## y grid distance in km
dny = floor(threshold/dky)  ## Number of grid boxes to search in y direction

## Create stencil library
stencils = list()
for (i in 1L:ny) {
  
  y1 = min(i + dny,ny)
  y0 = max(i - dny,1L)
  
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
                             cbind(lon1[grid[,1L]+dnx+1L],lat[grid[,2L]]))/1000
  
  stencil = grid[dists <= threshold,]
  stencil = cbind(stencil, cos(lat[stencil[,2L]]*pi/180))
  stencils[[i]] = stencil
  
} ## i


## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude" )
time.dim = ncdim_def("time", time.units, time,
                     unlim = TRUE, calendar = calendar, longname = "Time")


######################
## Temperature file ##
######################

## Attributes
global.attributes = ncatt_get(ncid, 0L)
atts = global.attributes[! names(global.attributes) %in% "history"]
make_missing_value = ncid$var[[vard]]$make_missing_value

## Define variables
precd = ncid$var[[vard]]$prec
if (precd == "int")
  precd = "integer"
d.var = ncvar_def(vard, ncid$var[[vard]]$units,
                  list(lon.dim,lat.dim,time.dim),
                  if (make_missing_value) ncid$var[[vard]]$missval else NULL,
                  longname = ncid$var[[vard]]$longname,
                  prec = precd, compression = 5L)

## Create netCDF file
ncod = nc_create(doutfile, list(d.var))

## Write standard names
ncatt_put(ncod, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(ncod, "latitude", "standard_name", "latitude", prec = "text")
ncatt_put(ncod, "time", "standard_name", "time", prec = "text")

## Write axes
ncatt_put(ncod, "longitude", "standard_name", "X", prec = "text")
ncatt_put(ncod, "latitude" , "standard_name", "Y", prec = "text")
ncatt_put(ncod, "time"     , "standard_name", "T", prec = "text")

## Write global attributes
for (att in names(atts))
  ncatt_put(ncod, 0L, att, atts[[att]])
ncatt_put(ncod, 0L, "threshold", threshold)

## Write history
history = paste0(format(Sys.time(), "%FT%XZ%z"), ": ", "./smooth_precip.R ",
                 threshold, " ", dinfile, " ", pinfile, " ", 
                 doutfile, " ", poutfile)
if ("history" %in% names(global.attributes))
  history = paste(history, global.attributes$history, sep = "\n")
ncatt_put(ncod, 0L, "history", history)


########################
## Precipitation file ##
########################

## Attributes
global.attributes = ncatt_get(ncip, 0L)
atts = global.attributes[! names(global.attributes) %in% "history"]
make_missing_value = ncip$var[[varp]]$make_missing_value

## Define variables
precd = ncip$var[[varp]]$prec
if (precd == "int")
  precd = "integer"
p.var = ncvar_def(varp, ncip$var[[varp]]$units,
                  list(lon.dim,lat.dim,time.dim),
                  if (make_missing_value) ncip$var[[varp]]$missval else NULL,
                  longname = ncip$var[[varp]]$longname,
                  prec = precd, compression = 5L)

## Create netCDF file
ncop = nc_create(poutfile, list(p.var))

## Write standard names
ncatt_put(ncop, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(ncop, "latitude", "standard_name", "latitude", prec = "text")
ncatt_put(ncop, "time", "standard_name", "time", prec = "text")

## Write axes
ncatt_put(ncop, "longitude", "standard_name", "X", prec = "text")
ncatt_put(ncop, "latitude" , "standard_name", "Y", prec = "text")
ncatt_put(ncop, "time"     , "standard_name", "T", prec = "text")

## Write global attributes
for (att in names(atts))
  ncatt_put(ncop, 0L, att, atts[[att]])
ncatt_put(ncop, 0L, "threshold", threshold)

## Write history
history = paste0(format(Sys.time(), "%FT%XZ%z"), ": ", "./smooth_precip.R ",
                 threshold, " ", dinfile, " ", pinfile, " ", 
                 doutfile, " ", poutfile)
if ("history" %in% names(global.attributes))
  history = paste(history, global.attributes$history, sep = "\n")
ncatt_put(ncop, 0L, "history", history)


###############
## Smoothing ##
###############

## Loop over times
for (t in 1L:nt) {
  
  ## Initialize storage
  doutput = matrix(0L, nx, ny)
  poutput = matrix(0L, nx, ny)
  
  ## Load data
  dinput = ncvar_get(ncid, vard, c(1L,1L,t), c(-1L,-1L,1L))
  pinput = ncvar_get(ncip, varp, c(1L,1L,t), c(-1L,-1L,1L))
  
  ## Loop over indices
  for (x in 1:nx) {
    for (y in 1:ny) {
      
      ## Use stencil to expand object
      stencil = stencils[[y]]
      stencil[,1L] = stencil[,1L] + x
      if (nx < nxt) {
        stencil = stencil[1L <= stencil[,1L] & stencil[,1L] <= nx,]
      } else {
        stencil[,1L] = (stencil[,1L] - 1L) %% nx + 1L
      }
      doutput[x,y] = weighted.mean(dinput[stencil[,1L:2L]], 
                                   stencil[,3L], na.rm = TRUE)
      poutput[x,y] = max(pinput[stencil[,1L:2L]], na.rm = TRUE)
      
      
    } ## y
  } ## x
  
  ## Write data
  ncvar_put(ncod, vard, doutput, c(1L,1L,t), c(-1L,-1L,1L))
  ncvar_put(ncop, varp, poutput, c(1L,1L,t), c(-1L,-1L,1L))
  
} ## t

## Close input files
nc_close(ncid)
nc_close(ncip)

## Close output file
nc_close(ncod)
nc_close(ncop)
