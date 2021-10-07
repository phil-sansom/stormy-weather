#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(ncdf4)

## Load source
source("src/image.R")
source("src/invertlat.R")
source("src/lonflip.R")

## Parse arguments
args = commandArgs(TRUE)
infile  = as.character(args[1])
outfile = as.character(args[2])

## Open file
nci = nc_open(infile)

## Get dimensions
lon0 = as.numeric(nci$dim$lon$vals)
lat0 = as.numeric(nci$dim$lat$vals)
time = ncvar_get(nci, "time")
time.units = ncatt_get(nci, "time", "units")$value
calendar = ncatt_get(nci, "time", "calendar")$value

nx = length(lon0)
ny = length(lat0)
nt = length(time)

## Transform dimensions, if necessary
fliplon = any(lon0[1] > 180)
fliplat = lat0[1] > lat0[2]
if (fliplon) {
  lon = lonflip(matrix(0, nx, ny), lon)$lon
} else {
  lon = lon0
}
if (fliplat) {
  lat = rev(lat0)
} else {
  lat = lat0
}

## Extract attributes
global.attributes = ncatt_get(nci, 0)

## Initialize storage
output = array(0, c(nx,ny,nt))

## Loop over times
for (t in 1:nt) {
  
  ## Load data
  input = ncvar_get(nci, start = c(1,1,t), count = c(-1,-1,1))

  ## Fix lon and lat if required  
  if (fliplon) {
    buffer = lonflip(input, lon0)$x
    input = buffer$x
  } else {
    lon = lon0
  }
  if (fliplat) {
    input = invertlat(input)
  }
  
  ## Duplicate data to cope with meridian crossing events
  buffer = rbind(input,input,input)
  
  ## Expand grid to take care of edge cases
  xx = matrix(0, nrow(buffer) + 2, ncol(buffer) + 2)
  xx[2:(nrow(buffer)+1),2:(ncol(buffer)+1)] = buffer
  
  ## Buffer dimensions
  nx1 = nrow(xx)

  ## Identify events
  nevents = 0
  events = list()
  points = which(xx > 0)
  npoints = length(points)
  while(npoints > 0) {
    
    ## Increment counter
    nevents = nevents + 1
    
    ## Starting point
    x0 = points[npoints] %% nx1
    y0 = points[npoints] %/% nx1 + 1
    
    ## Find all points connected to the starting point
    current = find.connected.component(xx, x0, y0)
    event = list()
    event$points = which(current > 0)
    event$x = event$points %% nx1 - 1
    event$y = event$points %/% nx1 #+ 1 - 1
    
    ## Store event
    events[[nevents]] = event
    
    ## Remove found points from master list
    mask = which(points %in% event$points)
    points = points[-mask]
    npoints = length(points)
    
  }
  
  ## Remove events not partially in central region
  i = 1
  while(i <= length(events)) {
    if (all(events[[i]]$x > 2*nx) | any(events[[i]]$x <= nx)) {
      events = events[-i]
    } else {
      events[[i]]$x = (events[[i]]$x - 1) %% nx + 1
      events[[i]]$points = (events[[i]]$y - 1)*nx + events[[i]]$x
      i = i + 1
    }
  }
  nevents = length(events)
  
  ## Store events
  result = matrix(0, nx, ny)
  for (i in 1:nevents)
    result[events[[i]]$points] = i
  output[,,t] = result
  
} ## t

## Close input file
nc_close(nci)


##################
## Write events ##
##################

## Attributes
atts = global.attributes[! names(global.attributes) %in% "history"]
make_missing_value = nci$var[[1]]$make_missing_value

## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude" )
time.dim = ncdim_def("time", time.units, time,
                     unlim = TRUE, calendar = calendar, longname = "Time")

## Define variables
events.var = ncvar_def("events", "", list(lon.dim,lat.dim,time.dim),
                      longname = "Extreme events", 
                      prec = "short",
                      compression = 5)

## Create netCDF file
nco = nc_create(outfile, list(events.var))

## Write data
ncvar_put(nco, "events", output)

## Write standard names
ncatt_put(nco, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(nco, "latitude" , "standard_name", "latitude" , prec = "text")
ncatt_put(nco, "time"     , "standard_name", "time"     , prec = "text")

## Write global attributes
for (att in names(atts))
  ncatt_put(nco, 0, att, atts[[att]])

## Write history
history = paste0(format(Sys.time(), "%FT%XZ%z"), ": ", "./label_events.R ",
                 infile, " ", outfile)
if ("history" %in% names(global.attributes))
  history = paste(history, global.attributes$history, sep = "\n")
ncatt_put(nco, 0, "history", history)

## Close output file
nc_close(nco)
