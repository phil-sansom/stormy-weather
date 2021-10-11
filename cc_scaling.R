#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(ncdf4)

## Load source
source("src/invertlat.R")
source("src/lonflip.R")

## Parse arguments
args = commandArgs(TRUE)
tfiles  = as.character(args[1])
pfiles  = as.character(args[2])
outfile = as.character(args[3])

## Read file lists
tlist = scan(tfiles, character(), -1, quiet = TRUE)
plist = scan(pfiles, character(), -1, quiet = TRUE)

## Quantile to compute
prob = 0.99

## Max memory to use: Defaults to 16GB
memory.to.use = 16*1024*1024*1024 

## Number of bins
nb = 12

## Temperature smoothing
smoothing = 0

## Threshold
threshold = 0.1/1000

## Read dimensions
nc  = nc_open(tlist[1])
lon0 = nc$dim$longitude$vals
lat0 = nc$dim$latitude$vals
calendar   = nc$dim$time$calendar
time.units = nc$dim$time$units
varname    = nc$var[[1]]$longname
units      = nc$var[[1]]$units
nc_close(nc)

nx = length(lon0)
ny = length(lat0)

## Determine if we need to transform dimensions
fliplon = any(lon0 >= 180)
fliplat = lat0[1] > lat0[2]

## Transform dimensions
lon = if (fliplon) lonflip(matrix(0, length(lon0), length(lat0)), lon0)$lon else lon0
lat = if (fliplat) rev(lat0) else lat0

## Read times
time = numeric()
for (file in tlist) {
  
  nc   = nc_open(file)
  time = c(time,nc$dim$time$vals)
  nc_close(nc)
  
} ## file
nt = length(time)

## Split into chunks
total.size = as.double(nx)*as.double(ny)*as.double(nt)*8*ifelse(smoothing > 0, 4, 3)
n.chunks = ceiling(total.size/memory.to.use)
chunk.size = ceiling(ny/n.chunks)
chunks = data.frame(
  start = seq(0, n.chunks - 1, 1)*chunk.size + 1,
  count = rep(chunk.size, n.chunks)
)
chunks$count[n.chunks] = ny - (n.chunks - 1)*chunk.size

## Initialize storage
temp   = array(NA, c(nx,ny,nb), 
               dimnames = list(longitude = lon, latitude = lat, bin = 1:nb))
precip = array(NA, c(nx,ny,nb), 
               dimnames = list(longitude = lon, latitude = lat, bin = 1:nb))

## Loop over chunks
for (i in 1:n.chunks) {
  
  start = chunks$start[i]
  count = chunks$count[i]
  temp0 = array(NA, c(nx,count,nt))
  precip0 = array(NA, c(nx,count,nt))
  
  ## Initialize time counter
  t1 = 1

  ## Loop over files
  for (j in 1:length(tlist)) {
    
    ## Print status
    print(paste("Chunk",i,"File",j))
    
    ## Open connections
    nct = nc_open(tlist[j])
    ncp = nc_open(plist[j])
    
    ## Get times
    ntj = nct$dim$time$len

    ## Load data
    mask = seq(t1, t1 + ntj - 1, 1)
    buffer = ncvar_get(nct, start = c(1,start,1), count = c(nx,count,ntj))
    temp0  [,,mask] = buffer
    buffer = ncvar_get(ncp, start = c(1,start,1), count = c(nx,count,ntj))
    precip0[,,mask] = buffer

    ## Close connections
    nc_close(nct)
    nc_close(ncp)
    rm(buffer)
    gc()

    ## Increment time counter
    t1 = t1 + ntj
  
  } ## j
  precip0[precip0 < 0] = 0
  gc()

  ## Smooth temp data
  if (smoothing > 0) {
    print(paste("Smoothing chunk", i))
    temps = temp0
    temp0[,,] = 0
    for (t in (1+smoothing):(nt-smoothing))
      for (s in -smoothing:+smoothing)
        temp0[,,t] = temp0[,,t] + temps[,,t + s]
    temp0 = temp0 / (2*smoothing + 1)
    rm(temps)
    gc()
  }

  ## Apply threshold
  precip0[precip0 < threshold] = NA
  gc()
  
  ## Bin data
  print(paste("Binning chunk", i))
  for (k in 1:nx) {
    for (l in 1:count) {
      mask = order(temp0[k,l,(1+smoothing):(nt-smoothing)]) + smoothing
      temp1 = temp0[k,l,mask]
      precip1 = precip0[k,l,mask]
      mask = !is.na(precip1)
      temp1 = temp1[mask]
      precip1 = precip1[mask]
      breaks = round(seq(0, length(precip1), length.out = nb + 1))
      for (m in 1:nb) {
        slice = (breaks[m]+1):breaks[m+1]
        temp[k,start+l-1,m] = mean(temp1[slice], na.rm = TRUE)
        precip[k,start+l-1,m] = quantile(precip1[slice], probs = prob, na.rm = TRUE)
      } ## m
    } ## l
  } ## k

  rm(precip0,temp0,temp1,precip1)
  gc()

} ## i

## Transform data
print("Transforming data...")
precip = 1000*precip
if (fliplon) {
  for (i in 1:nb) {
    temp[,,i] = lonflip(temp[,,i], lon0)$x
    precip[,,i] = lonflip(precip[,,i], lon0)$x
  } ## i
}
if (fliplat) {
  for (i in 1:nb) {
    temp[,,i] = invertlat(temp[,,i])
    precip[,,i] = invertlat(precip[,,i])
  } ## i
}


#######################
## Write climatology ##
#######################

print("Writing data...")

## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude")
bin.dim  = ncdim_def("bin", "number", 1:nb, longname = "Bin")
time.dim = ncdim_def("time", time.units, (time[2] + time[nt])/2,
                    unlim = TRUE, calendar = calendar, longname = "Time")
nv.dim = ncdim_def("bounds", "", 1:2, create_dimvar = FALSE)

## Define variables
temp.var  = ncvar_def("temp", "K", list(lon.dim,lat.dim,bin.dim,time.dim),  
                     longname = "Temperature", prec = "double",
                     compression = 5)
precip.var = ncvar_def("precip", "mm", list(lon.dim,lat.dim,bin.dim,time.dim),  
                     longname = "Precipitation", prec = "double",
                     compression = 5)
clim.var  = ncvar_def("climatology_bounds", "", list(nv.dim,time.dim),
                     prec = "double")


## Create netCDF file
nc = nc_create(outfile, list(clim.var, temp.var, precip.var))

## Write data
ncvar_put(nc, "temp", temp)
ncvar_put(nc, "precip", precip)
ncvar_put(nc, "climatology_bounds", c(time[1],time[nt]))

## Write description
ncatt_put(nc, 0, "Conventions", "CF-1.8", prec = "text")

## Write standard names
ncatt_put(nc, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(nc, "latitude" , "standard_name", "latitude" , prec = "text")
ncatt_put(nc, "time"     , "standard_name", "time"     , prec = "text")

## Write climatological time axis
ncatt_put(nc, "time", "climatology", "climatology_bounds", prec = "text")
ncatt_put(nc, "temp", "cell_methods", "time: mean", prec = "text")
ncatt_put(nc, "precip", "cell_methods", "time: quantile", prec = "text")

nc_close(nc)
