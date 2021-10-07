#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(ncdf4)

## Load source
source("src/invertlat.R")
source("src/lonflip.R")

## Parse arguments
args = commandArgs(TRUE)
nargs = length(args)
infiles = as.character(args[1:(nargs-1)])
outfile = as.character(args[nargs])

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
varname    = nc$var[[1]]$longname
units      = nc$var[[1]]$units
nc_close(nc)

np = length(probs)
nx = length(lon)
ny = length(lat)

## Read times
time = numeric()
for (file in infiles) {
  
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
  for (file in infiles) {
    
    ## Print status
    print(file)
    
    ## Open connection
    nc = nc_open(file)
    ntj = nc$dim$time$len

    ## Load data
    mask = seq(t1, t1 + ntj - 1, 1)
    chunk[,,mask] = ncvar_get(nc, start = c(1,start,1), count = c(nx,count,ntj))

    ## Close connection
    nc_close(nc)
    
    ## Increment time counter
    t1 = t1 + ntj
  
  } ## file
  
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

## Transform data, if necessary
if (any(lon > 180)) {
  means = lonflip(means, lon)$x
  sds   = lonflip(sds  , lon)$x
  for (p in 1:np)
    quantiles[,,p] = lonflip(quantiles[,,p], lon)$x
  lon = lonflip(matrix(0, nx, ny), lon)$lon
}
if (lat[1] > lat[2]) {
  means = invertlat(means)
  sds   = invertlat(sds)
  for (p in 1:np)
    quantiles[,,p] = invertlat(quantiles[,,p])
  lat = rev(lat)
}


#######################
## Write climatology ##
#######################

## Make climatology attributes
interval = time[2] - time[1]
climatology.time = median(c(time,time[length(time)] + interval))
climatology.bounds = c(time[1],time[length(time)] + interval)

## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude")
time.dim = ncdim_def("time", time.units, climatology.time,
                     unlim = TRUE, calendar = calendar, longname = "Time")
prob.dim = ncdim_def("probability", "", probs, longname = "Probability")
nv.dim   = ncdim_def("bounds", "", 1:2, create_dimvar = FALSE)

## Define variables
means.var = ncvar_def("mean", units, list(lon.dim,lat.dim,time.dim),  
                     longname = "Mean", prec = "double",
                     compression = 5)
sds.var   = ncvar_def("sd", units, list(lon.dim,lat.dim,time.dim),  
                     longname = "Standard deviation", prec = "double",
                     compression = 5)
quant.var = ncvar_def("quantiles", units, 
                     list(lon.dim,lat.dim,prob.dim,time.dim),  
                     longname = "Quantiles", prec = "double",
                     compression = 5)
clim.var  = ncvar_def("climatology_bounds", "", list(nv.dim,time.dim),
                     prec = "double", compression = 5)

## Create netCDF file
nc = nc_create(outfile, list(clim.var,means.var,sds.var,quant.var))

## Write data
ncvar_put(nc, "mean"     , means    )
ncvar_put(nc, "sd"       , sds      )
ncvar_put(nc, "quantiles", quantiles)
ncvar_put(nc, "climatology_bounds", climatology.bounds)

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
ncatt_put(nc, "sd"  , "cell_methods", "time: standard_deviation", prec = "text")
ncatt_put(nc, "quantiles", "cell_methods", "time: quantile", prec = "text")

nc_close(nc)
