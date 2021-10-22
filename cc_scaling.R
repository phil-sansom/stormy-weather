#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(optparse)
library(ncdf4)

## Load source
source("src/invertlat.R")
source("src/lonflip.R")

## Optional arguments
option_list = list(
  make_option(c("--quantile","-q"), action = "store", type = "double",
              help = "Quantile of precipitation to compute [default: 0.99]",
              default = 0.99),
  make_option(c("--smoothing","-s"), action = "store", type = "integer",
              help = "Number of time steps to smooth temperatures over [default: 0]",
              default = 0),
  make_option(c("--threshold","-t"), action = "store", type = "double",
              help = "Threshold for wet days (mm) [default: 0.1]",
              default = 0.1),
  make_option(c("--nbins","-n"), action = "store", type = "integer",
              help = "Number of equal size bins to use [default: 12]",
              default = 12),
  make_option(c("--compression","-c"), action = "store", type = "integer",
              help = "Compression level to use (0-9) [default: 5]",
              default = 5),
  make_option(c("--binsize","-b"), action = "store", type = "integer",
              help = "Minimum number of obs in each bin (overrides nbins)"),
  make_option(c("--memory","-m"), action = "store", type = "integer",
              help = "Maximum memory to use (in MB)")
)

## Argument parser
parser = OptionParser(usage = "Usage: %prog [OPTION]... TEMP_FILES PRECIP_FILES OUTFILE",
                      option_list = option_list,
                      description = "Clausius-Clapeyron scaling calculations.\n\nOperands:\n\tTEMP_FILES\n\t\tList of temperature files to read\n\tPRECIP_FILES\n\t\tList of precipitation files to read\n\tOUTFILE\n\t\tFile to write output to")

## Parser arguments
argv = parse_args(parser, positional_arguments = 3)
opts = argv$options
args = argv$args
smoothing = opts$smoothing
if (opts$compression == 0)
  opts$compression = NA

## Read file lists
tlist = scan(args[1], character(), -1, quiet = TRUE)
plist = scan(args[2], character(), -1, quiet = TRUE)
n.files = length(plist)

## Read dimensions
nc = nc_open(plist[1])
lon0 = nc$dim$longitude$vals
lat0 = nc$dim$latitude$vals
calendar   = nc$dim$time$calendar
time.units = nc$dim$time$units
precip.units = nc$var[[1]]$units
nc_close(nc)

nx = length(lon0)
ny = length(lat0)

## Determine if we need to transform dimensions
fliplon = any(lon0 >= 180)
fliplat = lat0[1] > lat0[2]

## Transform dimensions
if (fliplon) {
  lon = lonflip(matrix(0, length(lon0), length(lat0)), lon0)$lon 
} else {
  lon = lon0
}
lat = if (fliplat) rev(lat0) else lat0

## Read times
time = numeric()
for (file in plist) {
  
  nc   = nc_open(file)
  time = c(time,nc$dim$time$vals)
  nc_close(nc)
  
} ## file
nt = length(time)
dt = time[2] - time[1]
if (exists("binsize", opts)) {
  nb = floor((nt - 2*smoothing)/opts$binsize)
} else {
  nb = opts$nbins
}

## Split into chunks
if (exists("memory", opts)) {
  row.size   = nx*nt*8/1024/1024
  chunk.size = floor(opts$memory/row.size/2)
  n.chunks   = ceiling(ny/chunk.size)
} else {
  chunk.size = ny
  n.chunks = 1
}
chunks = data.frame(
  start = seq(0, n.chunks - 1, 1)*chunk.size + 1,
  count = rep(chunk.size, n.chunks)
)
chunks$count[n.chunks] = ny - (n.chunks - 1)*chunk.size
if (fliplat) {
  chunks$start = rev(chunks$start)
  chunks$count = rev(chunks$count)
}

## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude")
bin.dim  = ncdim_def("bin", "number", 1:nb, longname = "Bin")
time.dim = ncdim_def("time", time.units, (time[1] + time[nt] + dt)/2,
                     unlim = TRUE, calendar = calendar, longname = "Time")
nv.dim   = ncdim_def("bounds", "", 1:2, create_dimvar = FALSE)

## Define variables
temp.var    = ncvar_def("temp", "K", list(lon.dim,lat.dim,bin.dim,time.dim),
                        missval = 9.9692099683868690e+36,
                        longname = "Temperature", prec = "double",
                        compression = opts$compression)
precip.var  = ncvar_def("precip", "mm", list(lon.dim,lat.dim,bin.dim,time.dim),
                        missval = 9.9692099683868690e+36,
                        longname = "Precipitation", prec = "double",
                        compression = opts$compression)
binsize.var = ncvar_def("binsize", "", list(lon.dim,lat.dim,bin.dim,time.dim),
                        missval = -2147483647,
                        longname = "Bin size", prec = "integer",
                        compression = opts$compression)
clim.var    = ncvar_def("climatology_bounds", "", list(nv.dim,time.dim),
                        prec = "double")
vars = list(clim.var, temp.var, precip.var, binsize.var)
if (exists("nbins", opts)) {
  nbins.var = ncvar_def("nbins", "", list(lon.dim,lat.dim,time.dim),
                        longname = "Number of bins", prec = "integer",
                        compression = opts$compression)
  vars = c(vars,list(nbins.var))
}

## Create netCDF file
nco = nc_create(args[3], vars)

## Write description
ncatt_put(nco, 0, "Conventions", "CF-1.8", prec = "text")

## Write standard names
ncatt_put(nco, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(nco, "latitude" , "standard_name", "latitude" , prec = "text")
ncatt_put(nco, "time"     , "standard_name", "time"     , prec = "text")

## Write climatological time axis
ncvar_put(nco, "climatology_bounds", c(time[1],time[nt] + dt))
ncatt_put(nco, "time"  , "climatology" , "climatology_bounds", prec = "text")
ncatt_put(nco, "temp"  , "cell_methods", "time: mean"        , prec = "text")
ncatt_put(nco, "precip", "cell_methods", "time: quantile"    , prec = "text")

## Write parameters
ncatt_put(nco, 0, "quantile" , opts$quantile , prec = "double" )
ncatt_put(nco, 0, "smoothing", opts$smoothing, prec = "integer")
ncatt_put(nco, 0, "threshold", opts$threshold, prec = "double" )
if (exists("binsize", opts)) {
  ncatt_put(nco, 0, "binsize", opts$binsize, prec = "integer")
} else {
  ncatt_put(nco, 0, "nbins"  , opts$nbins  , prec = "integer")
}

## Transform threshold
if (! precip.units %in% c("mm","cm","m"))
  warning("Precip units not recognised (mm,cm,m), specify threshold in native units")
if (precip.units == "cm")
  opts$threshold = opts$threshold/10
if (precip.units == "m")
  opts$threshold = opts$threshold/1000

## Loop over chunks
for (i in 1:n.chunks) {
  
  ## Initialize storage
  start = chunks$start[i]
  count = chunks$count[i]
  temp0   = array(NA, c(nx,count,nt))
  precip0 = array(NA, c(nx,count,nt))
  temp    = array(NA, c(nx,count,nb))
  precip  = array(NA, c(nx,count,nb))
  binsize = array(NA, c(nx,count,nb))
  nbins   = array(NA, c(nx,count))
  
  ## Initialize time counter
  t1 = 1

  ## Loop over files
  for (j in 1:n.files) {
    
    ## Print status
    print(paste0("Reading chunk ",i," of ",n.chunks,
                 ", file ",j," of ", n.files))
    
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
    mask1 = which(buffer < 0)
    buffer[mask1] = 0
    mask1 = which(buffer < opts$threshold)
    buffer[mask1] = NA
    precip0[,,mask] = buffer

    ## Close connections
    nc_close(nct)
    nc_close(ncp)
    rm(buffer,mask,mask1)
    gc()

    ## Increment time counter
    t1 = t1 + ntj
  
  } ## j
  gc()

  # ## Smooth temp data
  if (smoothing > 0) {
    print(paste("Smoothing chunk",i,"of",n.chunks))
    smoothing = opts$smoothing
    nn = 2*smoothing + 1
    buffer = array(0, c(nx,count,nn))
    for (k in 1:nt) {
      kk = (k - 1) %% nn + 1
      buffer[,,kk] = 0
      mask = max(1,k - smoothing):min(k + smoothing,nt)
      for (l in mask)
        buffer[,,kk] = buffer[,,kk] + temp0[,,l]/length(mask)
      if (k > smoothing)
        temp0[,,k - smoothing] = buffer[,,(k - smoothing - 1) %% nn + 1]
    } ## k
    for (k in (nt - smoothing + 1):nt)
      temp0[,,k] = buffer[,,(k - 1) %% nn + 1]
    rm(mask,buffer,nn,kk)
    gc()
  }

  ## Bin data
  print(paste("Binning chunk",i,"of",n.chunks))
  for (k in 1:nx) {
    for (l in 1:count) {
      mask    = order(temp0[k,l,(1 + smoothing):(nt - smoothing)]) + smoothing
      temp1   = temp0  [k,l,mask]
      precip1 = precip0[k,l,mask]
      mask    = !is.na(precip1)
      temp1   = temp1  [mask]
      precip1 = precip1[mask]
      if (exists("binsize", opts)) {
        nbkl = floor(length(precip1)/opts$binsize)
      } else {
        nbkl = nb
      }
      nbins[k,l] = nbkl
      if (nbkl > 0) {
        breaks = round(seq(0, length(precip1), length.out = nbkl + 1))
        for (m in 1:nbkl) {
          slice = (breaks[m] + 1):breaks[m + 1]
          temp  [k,l,m] = mean(temp1[slice], na.rm = TRUE)
          precip[k,l,m] = quantile(precip1[slice], probs = opts$quantile, 
                                   na.rm = TRUE)
        } ## m
        binsize[k,l,1:nbkl] = diff(breaks)
      } ## nbkl > 0
    } ## l
  } ## k
  rm(precip0,temp0,temp1,precip1,mask,slice)
  gc()
  
  ## Transform data
  if (precip.units == "cm")
    precip = 10*precip
  if (precip.units == "m")
    precip = 1000*precip
  if (fliplon) {
    for (k in 1:nb) {
      temp   [,,k] = lonflip(temp   [,,k], lon0)$x
      precip [,,k] = lonflip(precip [,,k], lon0)$x
      binsize[,,k] = lonflip(binsize[,,k], lon0)$x
    } ## i
    nbins = lonflip(nbins, lon0)$x
  }
  if (fliplat & count > 1) {
    for (k in 1:nb) {
      temp   [,,k] = invertlat(temp   [,,k])
      precip [,,k] = invertlat(precip [,,k])
      binsize[,,k] = invertlat(binsize[,,k])
    } ## i
    nbins = invertlat(nbins)
  }
  
  ## Write data
  print(paste("Writing chunk",i,"of",n.chunks))
  if (fliplat) {
    ncvar_put(nco, "temp", temp, 
              start = c(1,ny - start - count + 2,1,1), count = c(nx,count,nb,1))
    ncvar_put(nco, "precip", precip, 
              start = c(1,ny - start - count + 2,1,1), count = c(nx,count,nb,1))
    ncvar_put(nco, "binsize", binsize, 
              start = c(1,ny - start - count + 2,1,1), count = c(nx,count,nb,1))
    if (exists("nbins", opts))
      ncvar_put(nco, "nbins", nbins, 
                start = c(1,ny - start - count + 2,1), count = c(nx,count,1))
  } else {
    ncvar_put(nco, "temp", temp, 
              start = c(1,start,1,1), count = c(nx,count,nb,1))
    ncvar_put(nco, "precip", precip, 
              start = c(1,start,1,1), count = c(nx,count,nb,1))
    ncvar_put(nco, "binsize", binsize, 
              start = c(1,start,1,1), count = c(nx,count,nb,1))
    if (exists("nbins", opts))
      ncvar_put(nco, "nbins", nbins, 
                start = c(1,start,1), count = c(nx,count,1))
  }
  
} ## i

## Close outout file
nc_close(nco)
