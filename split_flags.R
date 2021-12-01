#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(optparse)
library(ncdf4)

## Optional arguments
option_list = list(
  make_option(c("--compression","-c"), action = "store", type = "integer",
              help = "Compression level to use (0-9) [default: 5]",
              default = 5)
)

## Argument parser
parser = OptionParser(usage = "Usage: %prog [OPTION]... INFILE",
                      option_list = option_list,
                      description = "Split flag file into individual components.\n\nOperands:\n\tINFILE\n\t\tFile produced by flags.R")

## Parser arguments
argv = parse_args(parser, positional_arguments = 1)
opts = argv$options
args = argv$args
if (opts$compression == 0)
  opts$compression = NA

## Open files
nci = nc_open(args[1])

## Get dimensions
lon  = as.numeric(nci$dim$longitude$vals)
lat  = as.numeric(nci$dim$latitude$vals)
time = as.numeric(nci$dim$time$vals)
nx = length(lon)
ny = length(lat)
nt = length(time)
time.units = nci$dim$time$units
calendar   = nci$dim$time$calendar

## Get data
storm_type = ncvar_get(nci, "storm_type")
flag_values = ncatt_get(nci, "storm_type", "flag_values")$value
flag_meanings = ncatt_get(nci, "storm_type", "flag_meanings")$value
flag_meanings = strsplit(flag_meanings, " ")[[1]]

## Close input files
nc_close(nci)

## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude")
time.dim = ncdim_def("time", time.units, time,
                     unlim = TRUE, calendar = calendar, longname = "Time")

for (i in 1:length(flag_meanings)) {

  ## Output file
  flag = flag_meanings[i]
  flag_lwr = tolower(flag)
  if (!dir.exists(file.path(dirname(args[1]),flag_lwr)))
    dir.create(file.path(dirname(args[1]),flag_lwr))
  outfile = file.path(dirname(args[1]),flag_lwr,basename(args[1]))
  
  ## Define variables
  flag.var = ncvar_def(flag_lwr, "", list(lon.dim,lat.dim,time.dim),
                       missval = -127, longname = flag, prec = "byte",
                       compression = 5)
  
  ## Create netCDF file
  nco = nc_create(outfile, flag.var)
  
  ## Write description
  ncatt_put(nco, 0, "Conventions", "CF-1.9", prec = "text")
  
  ## Write standard names
  ncatt_put(nco, "longitude", "standard_name", "longitude", prec = "text")
  ncatt_put(nco, "latitude" , "standard_name", "latitude" , prec = "text")
  ncatt_put(nco, "time"     , "standard_name", "time"     , prec = "text")
  
  ## Write axes
  ncatt_put(nco, "longitude", "axis", "X", prec = "text")
  ncatt_put(nco, "latitude" , "axis", "Y", prec = "text")
  ncatt_put(nco, "time"     , "axis", "T", prec = "text")
  
  ## Write data
  ncvar_put(nco, flag_lwr, !is.na(storm_type) & storm_type == flag_values[i])
  
  ## Close output file
  nc_close(nco)
  
} ## i

