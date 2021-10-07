#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(ncdf4)

## Parse arguments
args = commandArgs(TRUE)
varid   = as.character(args[1])
infiles = as.character(args[2])
outfile = as.character(args[3])

## Read file list
file.list = scan(infiles, character(), -1, quiet = TRUE)

## Quantiles to compute
probs = c(0.005,0.01,0.02,0.025,0.05,0.10,0.25,0.50,
          0.75,0.90,0.95,0.975,0.98,0.99,0.995)

## Max memory to use: Defaults to 16GB
memory.to.use = 16*1024*1024*1024 

## Read dimensions
nc  = nc_open(file.list[1])
lon = nc$dim$longitude$vals
lat = nc$dim$latitude$vals
calendar   = nc$dim$time$calendar
time.units = nc$dim$time$units
varname    = nc$var[[varid]]$longname
units      = nc$var[[varid]]$units
nc_close(nc)

np = length(probs)
nx = length(lon)
ny = length(lat)

## Read times
time = numeric()
for (file in file.list) {
  
  nc   = nc_open(file)
  time = c(time,nc$dim$time$vals)
  nc_close(nc)
  
} ## file
nt = length(time)

## Split into chunks
total.size = as.double(nx)*as.double(ny)*as.double(nt)*8
n.chunks = ceiling(total.size/memory.to.use)
chunk.size = ceiling(ny/n.chunks)
chunks = data.frame(
  start = seq(0, n.chunks - 1, 1)*chunk.size + 1,
  count = rep(chunk.size, n.chunks)
)
chunks$count[n.chunks] = ny - (n.chunks - 1)*chunk.size

## Initialize storage
means     = array(NA, c(nx,ny), list(longitude = lon, latitude = lat))
sds       = array(NA, c(nx,ny), list(longitude = lon, latitude = lat))
quantiles = array(NA, c(nx,ny,np), list(longitude = lon, latitude = lat,
                                        quantile = probs))

## Loop over chunks
for (i in 1:n.chunks) {
  
  ## Print status
  print(paste("Chunk:", i))
  
  start = chunks$start[i]
  count = chunks$count[i]
  chunk = array(NA, c(nx,count,nt))
  
  ## Initialize time counter
  t1 = 1

  ## Loop over files
  for (j in 1:length(file.list)) {
    
    ## Print status
    print(file.list[j])
    
    ## Open connection
    nc = nc_open(file.list[j])
    ntj = nc$dim$time$len

    ## Load data
    mask = seq(t1, t1 + ntj - 1, 1)
    chunk[,,mask] = ncvar_get(nc, varid, c(1,start,1), c(nx,count,ntj))

    ## Close connection
    nc_close(nc)
    
    ## Increment time counter
    t1 = t1 + ntj
  
  } ## j
  
  ## Compute climatologies and store
  mask = seq(start, start + count - 1, 1)
  means    [,mask] = apply(chunk, c(1,2), mean, na.rm = TRUE)
  sds      [,mask] = apply(chunk, c(1,2), sd  , na.rm = TRUE)
  buffer = apply(chunk, c(1,2), quantile, probs = probs , na.rm = TRUE)
  if (np > 1) {
    for (j in 1:np)
      quantiles[,mask,j] = buffer[j,,]
  } else {
    quantiles[,mask,] = buffer
  }
  
} ## i

rm(chunk,buffer)


#######################
## Write climatology ##
#######################

## Define dimensions
lon_nc  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat_nc  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude")
time_nc = ncdim_def("time", time.units, (time[1] + time[2])/2,
                    unlim = TRUE, calendar = calendar, longname = "Time")
prob_nc = ncdim_def("probability", "", probs, longname = "Probability")
nv_nc = ncdim_def("bounds", "", 1:2, create_dimvar = FALSE)

## Define variables
means_nc = ncvar_def("mean", units, list(lon_nc,lat_nc,time_nc),  
                     longname = "Mean", prec = "float",
                     compression = 5)
sds_nc   = ncvar_def("sd", units, list(lon_nc,lat_nc,time_nc),  
                     longname = "Standard deviation", prec = "float",
                     compression = 5)
quant_nc = ncvar_def("quantiles", units, list(lon_nc,lat_nc,prob_nc,time_nc),  
                     longname = "Quantiles", prec = "float",
                     compression = 5)
clim_nc  = ncvar_def("climatology_bounds", "", list(nv_nc,time_nc),
                     prec = "double", compression = 5)


## Create netCDF file
nc = nc_create(outfile, list(clim_nc, means_nc, sds_nc, quant_nc))

## Write data
ncvar_put(nc, "mean"     , means    )
ncvar_put(nc, "sd"       , sds      )
ncvar_put(nc, "quantiles", quantiles)
ncvar_put(nc, "climatology_bounds", c(time[1],time[nt]))

## Write description
ncatt_put(nc, 0, "Conventions", "CF-1.8", prec = "text")
ncatt_put(nc, 0, "title", paste("Climatology:", varname), prec = "text")

## Write standard names
ncatt_put(nc, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(nc, "latitude" , "standard_name", "latitude" , prec = "text")
ncatt_put(nc, "time"     , "standard_name", "time"     , prec = "text")

## Write climatological time axis
ncatt_put(nc, "time", "climatology", "climatology_bounds", prec = "text")
ncatt_put(nc, "mean", "cell_methods", "time: mean", prec = "text")
ncatt_put(nc, "sd", "cell_methods", "time: standard_deviation", prec = "text")
ncatt_put(nc, "quantiles", "cell_methods", "time: quantile", prec = "text")

nc_close(nc)
