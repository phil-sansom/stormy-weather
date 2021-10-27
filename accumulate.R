#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(optparse)
library(ncdf4)

## Load source
source("src/history.R")

## Optional arguments
option_list = list(
  make_option("--start", action = "store", type = "integer",
              default = 1, help = "Start [default: 1]"),
  make_option("--step", action = "store", type = "integer",
              default = 6, help = "Time step [default: 6]"),
  make_option("--mask", action = "store", type = "character",
              default = "-2,-1,0,1,2,3", help = "Mask [default: -2,-1,0,1,2,3]"),
  make_option("--action", action = "store", type = "character",
              default = "max", help = "Action (max,sum,mean,min) [default: max]"),
  make_option("--previous", action = "store", type = "character",
              default = NA, help = "Previous input file"),
  make_option("--next", action = "store", type = "character",
              dest = "nextfile", default = NA, help = "Next input file"),
  make_option("--varname", action = "store", type = "character",
              default = NA, help = "Variable to read"),
  make_option("--compression", action = "store", type = "integer",
              default = 5, help = "Compression level to use (0-9) [default: 5]"),
  make_option("--pack", action = "store_true", type = "logical",
              default = FALSE, help = "Pack data")
)

## Argument parser
parser = OptionParser(usage = "Usage: %prog [OPTION]... INFILE OUTFILE",
                      option_list = option_list,
                      description = "Accumulation.\n\nOperands:\n\tINFILE\n\t\tInput file\n\tOUTFILE\n\t\tOutput file")

## Parser arguments
args = commandArgs(TRUE)
history = make.history("./accumulate.R", args)
argv = parse_args(parser, args = args, positional_arguments = 2)
opts = argv$options
args = argv$args
opts$mask = as.numeric(strsplit(opts$mask, ",")[[1]])
if (opts$compression == 0)
  opts$compression = NA

infile  = args[1]
outfile = args[2]

## Open connections
nci = nc_open(infile)
if (is.na(opts$varname))
  opts$varname = names(nci$var)[[1]]
if (!is.na(opts$previous)) {
  ncp = nc_open(opts$previous)
  ntp = ncp$dim$time$len
}
if (!is.na(opts$nextfile)) {
  ncn = nc_open(opts$nextfile)
  ntn = ncn$dim$time$len
}

## Read dimensions
time0 = nci$dim$time$vals

nx  = nci$dim$longitude$len
ny  = nci$dim$latitude$len
nt0 = nci$dim$time$len
nm  = length(opts$mask)

time.mask = seq(opts$start, nt0, opts$step)
time = time0[time.mask]
nt = length(time.mask)
if (time.mask[1] + opts$mask[1] < 1 & is.na(opts$previous))
  warning("Specify preceding file (--previous) or first time steps will not be computed correctly")
if (nt0 < time.mask[nt] + opts$mask[nm] & is.na(opts$nextfile))
  warning("Specify succeeding file (--next) or final time steps will not be computed correctly")
if (2*nt0 < nm)
  stop("Length of mask exceeds twice the number of time steps")

## Initialize storage
output = array(NA, c(nx,ny,nt))

## Loop over times
mask.min = min(opts$mask)
mask.max = max(opts$mask)
nm = mask.max - mask.min + 1
for (t in 1:nt) {

  count = nm
  countp = 0
  countn = 0
  buffer = array(NA, c(nx,ny,count))
  tt = time.mask[t]
  start = tt + mask.min
  if (start <= 0) {
    countp = 1 - start
    if (!is.na(opts$previous)) {
      startp = ntp + start - 1
      buffer[,,1:countp] = ncvar_get(ncp, opts$varname, start = c(1,1,startp), 
                                     count = c(nx,ny,countp))
    }
    start = 1
    count = nm - countp
  }
  if (nt0 < start + count - 1) {
    startn = 1
    countn = start + count - nt0 - 1
    if (!is.na(opts$nextfile)) {
      buffer[,,(nm - countn + 1):nm] = ncvar_get(ncp, opts$varname, 
                                                 start = c(1,1,startn), 
                                                 count = c(nx,ny,countn))
    }
    count = count - countn
  }
  buffer[,,(countp + 1):(nm - countn)] = 
    ncvar_get(nci, opts$varname, start = c(1,1,start), count = c(nx,ny,count))
  output[,,t] = switch(opts$action,
                       min  = apply(buffer, c(1,2), min , na.rm = TRUE),
                       max  = apply(buffer, c(1,2), max , na.rm = TRUE),
                       mean = apply(buffer, c(1,2), mean, na.rm = TRUE),
                       sum  = apply(buffer, c(1,2), sum , na.rm = TRUE))
  
} ## t

## Define dimensions
lon.dim  = ncdim_def("longitude", nci$dim$longitude$units, nci$dim$longitude$vals)
lat.dim  = ncdim_def("latitude" , nci$dim$latitude$units, nci$dim$latitude$vals)
time.dim = ncdim_def("time", nci$dim$time$units, time, unlim = TRUE, 
                     calendar = nci$dim$time$calendar)

## Define variables
if (opts$pack) {
  data.var = ncvar_def(opts$varname, nci$var[[opts$varname]]$units, 
                       list(lon.dim,lat.dim,time.dim),
                       missval = -32767, prec = "short")
} else {
  data.var = ncvar_def(opts$varname, nci$var[[opts$varname]]$units, 
                       list(lon.dim,lat.dim,time.dim),
                       missval = 9.9692099683868690e+36, prec = "double",
                       compression = opts$compression)
}
  
## Create netCDF file
nco = nc_create(outfile, data.var)

## Write attributes
global.attributes = ncatt_get(nci, 0)
if (length(global.attributes) > 0)
  for (i in 1:length(global.attributes))
    ncatt_put(nco, 0, names(global.attributes)[i], global.attributes[[i]])
lon.attributes = ncatt_get(nci, "longitude")
if (length(lon.attributes) > 0)
  for (i in 1:length(lon.attributes))
    ncatt_put(nco, "longitude", names(lon.attributes)[i], lon.attributes[[i]])
lat.attributes = ncatt_get(nci, "latitude")
if (length(lat.attributes) > 0)
  for (i in 1:length(lat.attributes))
    ncatt_put(nco, "latitude", names(lat.attributes)[i], lat.attributes[[i]])
time.attributes = ncatt_get(nci, "time")
if (length(time.attributes) > 0)
  for (i in 1:length(time.attributes))
    ncatt_put(nco, "time", names(time.attributes)[i], time.attributes[[i]])
var.attributes = ncatt_get(nci, opts$varname)
var.attributes = var.attributes[!names(var.attributes) %in% 
                                  c("_FillValue","_NoFill","missing_value",
                                    "valid_min","valid_max","valid_range",
                                    "scale_factor","add_offset")]
if (length(var.attributes) > 0)
  for (i in 1:length(var.attributes))
    ncatt_put(nco, opts$varname, names(var.attributes)[i], var.attributes[[i]])
ncatt_put(nco, opts$varname, "missing_value", data.var$missval)

## Write data
if (opts$pack) {
  n = 16
  scale.factor = (max(output) - min(output))/(2^n - 3)
  add.offset = 
    ((2^(n - 1) - 1)*(max(output) - min(output)) - max(output))/(2^n - 3)
  packed = round((output - add.offset)/scale.factor)
  ncvar_put(nco, opts$varname, packed)
  ncatt_put(nco, opts$varname, "scale_factor", scale.factor, prec = "double")
  ncatt_put(nco, opts$varname, "add_offset"  , add.offset  , prec = "double")
} else {
  ncvar_put(nco, opts$varname, output)
}

## Close connections
nc_close(nci)
nc_close(nco)
if (!is.na(opts$previous))
  nc_close(ncp)
if (!is.na(opts$nextfile))
  nc_close(ncn)
