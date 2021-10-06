#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(ncdf4)

## Parse arguments
args = commandArgs(TRUE)
varid   = as.character(args[1])
infile  = as.character(args[2])
outfile = as.character(args[3])

## Open connection
nci = nc_open(infile)

## Get dimensions
longitude = as.numeric(nci$dim$longitude$vals)
latitude  = as.numeric(nci$dim$latitude$vals)
time      = as.numeric(nci$dim$time$vals)

## Get attributes
attributes = list()
attributes$global = ncatt_get(nci, 0)
for (dim in names(nci$dim))
  attributes[[dim]] = ncatt_get(nci, dim)
for (var in names(nci$var))
  attributes[[var]] = ncatt_get(nci, var)

## Load data
data = ncvar_get(nci, varid)

## Close connection
nc_close(nci)

## Difference data
mask1 = seq(1, length(time) - 1, 2)
mask2 = seq(2, length(time), 2)
data[,,mask2] = data[,,mask2] - data[,,mask1]

## Define dimensions
lon.dim  = ncdim_def("longitude", attributes$longitude$units, longitude, 
                     longname = attributes$longitude$long_name)
lat.dim  = ncdim_def("latitude", attributes$latitude$units, latitude, 
                     longname = attributes$latitude$long_name)
time.dim = ncdim_def("time", attributes$time$units, time,
                     calendar = attributes$time$calendar,
                     longname = attributes$time$long_name)

## Define variables
data.var  = ncvar_def(varid, attributes[[varid]]$units, 
                      list(lon.dim,lat.dim,time.dim),  
                      longname = attributes[[varid]]$long_name,
                      compression = 5)

## Create netCDF file
nco = nc_create(outfile, data.var)

## Write global attributes
for (attribute in names(attributes$global))
  ncatt_put(nco, 0, attribute, attributes$global[[attribute]])

## Write data
ncvar_put(nco, varid, data)

## Close connection
nc_close(nco)
