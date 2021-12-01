#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(optparse)
library(ncdf4)

## Optional arguments
option_list = list(
  make_option(c("--none","-n"), action = "store_true", type = "logical",
              help = "Include instances with no storms", default = FALSE),
  make_option(c("--compression","-c"), action = "store", type = "integer",
              help = "Compression level to use (0-9) [default: 5]",
              default = 5)
)

## Argument parser
parser = OptionParser(usage = "Usage: %prog [OPTION]... CYCLONES FRONTS THUNDER OUTFILE",
                      option_list = option_list,
                      description = "Combine cyclones, fronts and thunderstorms for classification.\n\nOperands:\n\tCYCLONES\n\t\tFile containing cyclones\n\tFRONTS\n\t\tFile containing fronts\n\tTHUNDER\n\t\tFile containing thunderstorms\n\tOUTFILE\n\t\tFile to write output to")

## Parser arguments
argv = parse_args(parser, positional_arguments = 4)
opts = argv$options
args = argv$args
if (opts$compression == 0)
  opts$compression = NA

## Open files
ncc = nc_open(args[1])
ncf = nc_open(args[2])
nct = nc_open(args[3])

## Get dimensions
lon   = as.numeric(ncc$dim$longitude$vals)
lat   = as.numeric(ncc$dim$latitude$vals)
time  = as.numeric(ncc$dim$time$vals)
lonf  = as.numeric(ncf$dim$longitude$vals)
latf  = as.numeric(ncf$dim$latitude$vals)
timef = as.numeric(ncf$dim$time$vals)
lont  = as.numeric(nct$dim$longitude$vals)
latt  = as.numeric(nct$dim$latitude$vals)
timet = as.numeric(nct$dim$time$vals)
nx = length(lon)
ny = length(lat)
nt = length(time)
time.units = ncc$dim$time$units
calendar = ncc$dim$time$calendar

## Check dimensions
if (any(lonf != lon) | any(lont != lon))
  stop("Longitude dimensions don't match")
if (any(latf != lat) | any(latt != lat))
  stop("Latitude dimensions don't match")
if (any(timef != time) | any(timet != time))
  stop("Time dimensions don't match")

## Get data
cyclones = ncvar_get(ncc)
fronts   = ncvar_get(ncf)
thunder  = ncvar_get(nct)
cyclones = cyclones > 0
cyclones[is.na(cyclones)] = FALSE
fronts = fronts > 0
fronts[is.na(fronts)] = FALSE
thunder = thunder > 0
thunder[is.na(thunder)] = FALSE

## Close input files
nc_close(ncc)
nc_close(ncf)
nc_close(nct)

## Combine data
data = array(NA, c(nx,ny,nt))
if (opts$none) {
  data[!(cyclones | fronts | thunder)] = 0
} else {
  data[!(cyclones | fronts | thunder)] = NA
}
data[cyclones & !(fronts | thunder)] = 1
data[fronts & !(cyclones | thunder)] = 2
data[thunder & !(cyclones | fronts)] = 3
data[(cyclones & fronts) & !thunder] = 4
data[(cyclones & thunder) & !fronts] = 5
data[(fronts & thunder) & !cyclones] = 6
data[fronts & thunder & cyclones]    = 7

## Data dictionary
if (opts$none) {
  flag_values = c(0,1,2,3,4,5,6,7)
  flag_meanings = "N CO FO TO CF CT FT CFT"
} else {
  flag_values = c(1,2,3,4,5,6,7)
  flag_meanings = "CO FO TO CF CT FT CFT"
}

## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude")
time.dim = ncdim_def("time", time.units, time,
                     unlim = TRUE, calendar = calendar, longname = "Time")

## Define variables
storms.var = ncvar_def("storm_type", "", list(lon.dim,lat.dim,time.dim),
                       missval = -127, longname = "Storm type", prec = "byte",
                       compression = 5)

## Create netCDF file
nco = nc_create(args[4], storms.var)

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

## Write flags
ncatt_put(nco, "storm_type", "flag_values", flag_values, prec = "byte")
ncatt_put(nco, "storm_type", "flag_meanings", flag_meanings, prec = "text")

## Write data
ncvar_put(nco, "storm_type", data)

## Close output file
nc_close(nco)
