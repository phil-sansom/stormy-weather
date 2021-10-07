#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(ncdf4)

## Load source
source("src/invertlat.R")
source("src/lonflip.R")

## Parse arguments
args = commandArgs(TRUE)
thresh  = as.character(args[1])  ## Threshold file
infile  = as.character(args[2])  ## Input file
outfile = as.character(args[3])  ## Output file

## Read threshold file
nc = nc_open(thresh)
lont = as.numeric(nc$dim$lon$vals)
latt = as.numeric(nc$dim$lat$vals)
threshold = ncvar_get(nc)
nc_close(nc)

## Transform threshold file, if necessary
if (any(lont > 180)) {
  buffer = lonflip(threshold, lont)
  threshold = buffer$x
  lont = buffer$lon
}
if (latt[1] > latt[2]) {
  threshold = invertlat(threshold)
  latt = rev(latt)
}

## Open input file
nci = nc_open(infile)

## Get dimensions
lon0 = as.numeric(nci$dim$lon$vals)
lat0 = as.numeric(nci$dim$lat$vals)
times = ncvar_get(nci, "time")
time.units = ncatt_get(nci, "time", "units")$value
calendar = ncatt_get(nci, "time", "calendar")$value

nx = length(lon0)
ny = length(lat0)
nt = length(times)

## Transform dimensions, if necessary
fliplon = any(lon0 > 180)
fliplat = lat0[1] > lat0[2]
if (fliplon) {
  lon = lonflip(matrix(0, nx, ny), lon0)$lon
} else {
  lon = lon0
}
if (fliplat) {
  lat = rev(lat0)
} else {
  lat = lat0
}

## Check compatibility of threshold and data
if(!(identical(latt,lat) & identical(lont,lon)))
  stop("Dimensions of threshold and data do not match")

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
    input = lonflip(input, lon0)$x
  }
  if (fliplat) {
    input = invertlat(input)
  }
  
  output[,,t] = 1*(input >= threshold)
  
} ## t

## Close input file
nc_close(nci)


##################
## Write events ##
##################

## Attributes
atts = global.attributes[! names(global.attributes) %in% "history"]

## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude" )
time.dim = ncdim_def("time", time.units, times,
                     unlim = TRUE, calendar = calendar, longname = "Time")

## Define variables
events.var = ncvar_def("events", "", list(lon.dim,lat.dim,time.dim),
                      longname = "Extreme events", 
                      prec = "byte",
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
history = paste0(format(Sys.time(), "%FT%XZ%z"), ": ", "./extract_events.R ",
                 thresh, " ", infile, " ", outfile)
if ("history" %in% names(global.attributes))
  history = paste(history, global.attributes$history, sep = "\n")
ncatt_put(nco, 0, "history", history)

## Close output file
nc_close(nco)
