#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(ncdf4)

## Load source
source("src/invertlat.R")
source("src/lonflip.R")

## Parse arguments
args = commandArgs(TRUE)
maskid  = as.character(args[1])  ## Threshold variable
thresh  = as.character(args[2])  ## Threshold file
varid   = as.character(args[3])  ## Variable name
infile  = as.character(args[4])  ## Input file
outfile = as.character(args[5])  ## Output file

## Read threshold file
nc = nc_open(thresh)
threshold = ncvar_get(nc, maskid)
nc_close(nc)

## Open input file
nci = nc_open(infile)

## Get longitudes and latitudes
lon0 = as.numeric(nci$dim$lon$vals)
lat0 = as.numeric(nci$dim$lat$vals)
fliplon = lon0[1] < 0
fliplat = lat0[1] > lat0[2]
if (fliplon) {
  buffer = lonflip(threshold, lon)
  threshold = buffer$x
  lon = buffer$lon
} else {
  lon = lon0
}
if (fliplat) {
  threshold = invertlat(threshold)
  lat = rev(lat0)
} else {
  lat = lat0
}

## Get times
times = ncvar_get(nci, "time")
time_units = ncatt_get(nci, "time", "units")$value
calendar = ncatt_get(nci, "time", "calendar")$value

## Dimensions
nlon = length(lon)
nlat = length(lat)
nt = length(times)

## Extract attributes
global.attributes = ncatt_get(nci, 0)

## Initialize storage
output = array(0, c(nlon,nlat,nt))

## Loop over times
for (t in 1:nt) {
  
  ## Load data
  input = ncvar_get(nci, varid, c(1,1,t), c(-1,-1,1))

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
lon_nc  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat_nc  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude" )
time_nc = ncdim_def("time", time_units, times,
                    unlim = TRUE, calendar = calendar, longname = "Time")

## Define variables
events_nc = ncvar_def("events", "", list(lon_nc,lat_nc,time_nc), 0,
                      longname = "Extreme events", 
                      prec = "byte",
                      compression = 5)

## Create netCDF file
nco = nc_create(outfile, list(events_nc))

## Write data
ncvar_put(nco, "events", output)

## Write standard names
ncatt_put(nco, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(nco, "latitude", "standard_name", "latitude", prec = "text")
ncatt_put(nco, "time", "standard_name", "time", prec = "text")

## Write global attributes
for (att in names(atts))
  ncatt_put(nco, 0, att, atts[[att]])

## Write history
history = paste0(format(Sys.time(), "%FT%XZ%z"), ": ", "./extract_events.R ",
                 maskid, " ", thresh, " ", varid, " ", infile, " ", outfile)
if ("history" %in% names(global.attributes))
  history = paste(history, global.attributes$history, sep = "\n")
ncatt_put(nco, 0, "history", history)

## Close output file
nc_close(nco)
