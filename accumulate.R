#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(optparse)
library(ncdf4)

## Optional arguments
option_list = list(
  make_option("--start", action = "store", type = "integer",
              help = "Start [default: 1]",
              default = 0.99),
  make_option("--step", action = "store", type = "integer",
              help = "Time step [default: 6]",
              default = 0),
  make_option("--mask", action = "store", type = "character",
              help = "Mask [default: -2,-1,0,1,2,3]",
              default = "-2,-1,0,1,2,3"),
  make_option("--action", action = "store", type = "character",
              help = "Action (max,sum,mean,min) [default: max]",
              default = "max"),
  make_option("--previous", action = "store", type = "character",
              help = "Previous input file"),
  make_option("--next", action = "store", type = "character",
              dest = "nextfile", help = "Next input file"),
  make_option("--compression", action = "store", type = "integer",
              help = "Compression level to use (0-9) [default: 5]",
              default = 5),
  make_option("--pack", action = "store_true", type = "logical",
              help = "Pack data",
              default = FALSE)
)

## Argument parser
parser = OptionParser(usage = "Usage: %prog [OPTION]... INFILE OUTFILE",
                      option_list = option_list,
                      description = "Accumulation.\n\nOperands:\n\tINFILE\n\t\tInput file\n\tOUTFILE\n\t\tOutput file")

## Parser arguments
argv = parse_args(parser, positional_arguments = 2)
opts = argv$options
args = argv$args
opts$mask = as.numeric(strsplit(opts$mask, ",")[[1]])
if (opts$compression == 0)
  opts$compression = NA
infile  = args[1]
outfile = args[2]

## Open connections
nci = nc_open(infile)
if (exists("previous", opts)) {
  ncp = nc_open(opts$previous)
  ntp = ncp$dim$time$len
}
if (exists("nextfile", opts)) {
  ncn = nc_open(opts$nextfile)
  ntn = ncn$dim$time$len
}

## Read dimensions
lon  = nci$dim$longitude$vals
lat  = nci$dim$latitude$vals
time0 = nci$dim$time$vals
time.calendar = nci$dim$time$calendar
time.units    = nci$dim$time$units
var.name      = nci$var[[1]]$name
var.units     = nci$var[[1]]$units
var.longname  = nci$var[[1]]$longname

nx  = length(lon)
ny  = length(lat)
nt0 = length(time0)
nm  = length(opts$mask)

time.mask = seq(opts$start, nt0, opts$step)
time = time0[time.mask]
nt = length(time.mask)
if (time.mask[1] + opts$mask[1] < 1 & !exists("previous", opts))
  warning("Specify preceding file (--previous) or first time steps will not be computed")
if (nt < time.mask[nt] + opts$mask[nm] & !exists("nextfile", opts))
  warning("Specify succeeding file (--next) or final time steps will not be computed")
if (2*nt < nm)
  stop("Length of mask exceeds twice the number of time steps")

## Initialize storage
output = array(NA, c(nx,ny,nt))

## Loop over times
mask.min = min(opts$mask)
mask.max = max(opts$mask)
mask.range = mask.max - mask.min + 1
nm = length(opts$mask)
for (t in 1:nt) {
  print(t)
  
  count = nm
  countp = 0
  countn = 0
  buffer = array(NA, c(nx,ny,count))
  tt = time.mask[t]
  start = tt + mask.min
  if (start <= 0) {
    countp  = 1 - start
    if (exists("previous", opts)) {
      startp  = ntp + start - 1
      buffer[,,1:countp] = 
        ncvar_get(ncp, var.name, start = c(1,1,startp), count = c(nx,ny,countp))
    }
    start = 1
    count = nm - countp
  }
  if (nt0 < tt + count) {
    startn = 1
    countn = tt + count - nt0 - 1
    if (exists("nextfile", opts)) {
      buffer[,,(nm - countn + 1):nm] = 
        ncvar_get(ncp, var.name, start = c(1,1,startn), count = c(nx,ny,countn))
    }
    count = count - countn
  }
  buffer[,,(countp+1):(nm - countn)] = 
    ncvar_get(nci, var.name, start = c(1,1,start), count = c(nx,ny,count))
  output[,,t] = switch(opts$action,
                       min  = apply(buffer, c(1,2), min , na.rm = TRUE),
                       max  = apply(buffer, c(1,2), max , na.rm = TRUE),
                       mean = apply(buffer, c(1,2), mean, na.rm = TRUE),
                       sum  = apply(buffer, c(1,2), sum , na.rm = TRUE))
  
} ## t

## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude" )
time.dim = ncdim_def("time", time.units, time,
                     unlim = TRUE, calendar = time.calendar, longname = "Time")

## Define variables
if (opts$pack) {
  data.var = ncvar_def(var.name, var.units, list(lon.dim,lat.dim,time.dim),
                       missval = -32767, 
                       longname = var.longname, prec = "short")
} else {
  data.var = ncvar_def(var.name, var.units, list(lon.dim,lat.dim,time.dim),
                       missval = 9.9692099683868690e+36,
                       longname = var.longname, prec = "double",
                       compression = opts$compression)
}
  
## Create netCDF file
nco = nc_create(outfile, data.var)

## Write description
ncatt_put(nco, 0, "Conventions", "CF-1.8", prec = "text")

## Write standard names
ncatt_put(nco, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(nco, "latitude" , "standard_name", "latitude" , prec = "text")
ncatt_put(nco, "time"     , "standard_name", "time"     , prec = "text")

## Write data
if (opts$pack) {
  n = 16
  scale.factor = (max(output) - min(output))/(2^n - 3)
  add.offset = 
    ((2^(n - 1) - 1)*(max(output) - min(output)) - max(output))/(2^n - 3)
  packed = round((output - add.offset)/scale.factor)
  ncvar_put(nco, var.name, packed)
  ncatt_put(nco, var.name, "scale_factor", scale.factor, prec = "double")
  ncatt_put(nco, var.name, "add_offset"  , add.offset  , prec = "double")
} else {
  ncvar_put(nco, var.name, output)
}

## Close connections
nc_close(nci)
nc_close(nco)
if (exists("previous", opts))
  nc_close(ncp)
if (exists("next", opts))
  nc_close(ncn)
