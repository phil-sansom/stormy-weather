#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(optparse)
library(ncdf4)
suppressPackageStartupMessages(library(quantreg))

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
  make_option(c("--level","-l"), action = "store", type = "double",
              help = "Nominal confidence level required [default: 0.95]",
              default = 0.95),
  make_option(c("--compression","-c"), action = "store", type = "integer",
              help = "Compression level to use (0-9) [default: 5]",
              default = 5),
  make_option(c("--memory","-m"), action = "store", type = "integer",
              help = "Maximum memory to use (in MB)"),
  make_option("--mask", action = "store", type = "character",
              help = "List of mask files to apply")
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
if (length(tlist) != length(plist))
  stop("Lengths of file lists differ")
if (exists("mask", opts)) {
  mlist = scan(opts$mask, character(), -1, quiet = TRUE)
  if (length(tlist) != length(mlist))
    stop("Length of mask list differs from temp and precip lists")
}
n.files = length(tlist)

## Check dimensions
time = numeric()
for (i in 1:n.files) {
  
  nct = nc_open(tlist[i])
  lont = nct$dim$longitude$vals
  latt = nct$dim$latitude$vals
  timet = nct$dim$time$vals
  nc_close(nct)
  ncp = nc_open(plist[i])
  lonp = ncp$dim$longitude$vals
  latp = ncp$dim$latitude$vals
  timep = ncp$dim$time$vals
  nc_close(ncp)
  if (any(lont != lonp))
    stop(paste("Longitude dimensions don't match in file", i))
  if (any(latt != latp))
    stop(paste("Latitude dimensions don't match in file", i))
  if (any(timet != timep))
    stop(paste("Time dimensions don't match in file", i))
  if (i > 1) {
    if (any(lont != lont0))
      stop(paste("Longitude dimensions don't match between temp files", i-1, "and", i))
    if (any(latt != latt0))
      stop(paste("Latitude dimensions don't match between temp files",  i-1, "and", i))
    if (any(lonp != lonp0))
      stop(paste("Longitude dimensions don't match between precip files", i-1, "and", i))
    if (any(latp != latp0))
      stop(paste("Latitude dimensions don't match between precip files",  i-1, "and", i))
  }
  lont0 = lont; latt0 = latt; timet0 = timet
  lonp0 = lonp; latp0 = latp; timep0 = timep
  
  if (exists("mask", opts)) {
    ncm = nc_open(mlist[i])
    lonm = ncm$dim$longitude$vals
    latm = ncm$dim$latitude$vals
    timem = ncm$dim$time$vals
    nc_close(ncm)
    if (any(lont != lonm))
      stop(paste("Longitude dimensions don't match in mask file", i))
    if (any(latt != latm))
      stop(paste("Latitude dimensions don't match in mask file", i))
    if (any(timet != timem))
      stop(paste("Time dimensions don't match in mask file", i))
    if (i > 1) {
      if (any(lonm != lonm0))
        stop(paste("Longitude dimensions don't match between mask files", i-1, "and", i))
      if (any(latm != latm0))
        stop(paste("Latitude dimensions don't match between mask files",  i-1, "and", i))
    }
    lonm0 = lonm; latm0 = latm; timem0 = timem
  }
  
  time = c(time,timet)
  
} ## i
nt = length(time)
dt = time[2] - time[1]

## Read dimensions
nc = nc_open(tlist[1])
lon0 = nc$dim$longitude$vals
lat0 = nc$dim$latitude$vals
calendar   = nc$dim$time$calendar
time.units = nc$dim$time$units
temp.name  = nc$var[[1]]$name
temp.units = nc$var[[1]]$units
temp.longname = nc$var[[1]]$longname
nc_close(nc)

nc = nc_open(plist[1])
precip.name  = nc$var[[1]]$name
precip.units = nc$var[[1]]$units
precip.longname = nc$var[[1]]$longname
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

## Split into chunks
if (exists("mask", opts)) {
  nd = 3
} else {
  nd = 2
}
if (exists("memory", opts)) {
  row.size   = nx*nt*8/1024/1024
  chunk.size = floor(opts$memory/row.size/nd)
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

## Define variables
beta.var = ncvar_def("beta", "%/K", list(lon.dim,lat.dim),
                     missval = 9.9692099683868690e+36,
                     longname = "Fitted value", prec = "double",
                     compression = opts$compression)
lwr.var = ncvar_def("lower", "%/K", list(lon.dim,lat.dim),
                     missval = 9.9692099683868690e+36,
                     longname = "Lower bound", prec = "double",
                     compression = opts$compression)
upr.var = ncvar_def("upper", "%/K", list(lon.dim,lat.dim),
                     missval = 9.9692099683868690e+36,
                     longname = "Upper bound", prec = "double",
                     compression = opts$compression)
vars = list(beta.var, lwr.var, upr.var)

## Create netCDF file
nco = nc_create(args[3], vars)

## Write description
ncatt_put(nco, 0, "Conventions", "CF-1.9", prec = "text")

## Write standard names
ncatt_put(nco, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(nco, "latitude" , "standard_name", "latitude" , prec = "text")

## Write axes
ncatt_put(nco, "longitude", "axis", "X", prec = "text")
ncatt_put(nco, "latitude" , "axis", "Y", prec = "text")

## Write parameters
ncatt_put(nco, 0, "quantile" , opts$quantile , prec = "double" )
ncatt_put(nco, 0, "smoothing", opts$smoothing, prec = "integer")
ncatt_put(nco, 0, "threshold", opts$threshold, prec = "double" )
ncatt_put(nco, 0, "level"    , opts$level    , prec = "double" )

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
  if (exists("mask", opts))
    mask0   = array(NA, c(nx,count,nt))
  beta  = array(NA, c(nx,count))
  lower = array(NA, c(nx,count))
  upper = array(NA, c(nx,count))

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
    if (exists("mask", opts))
      ncm = nc_open(mlist[j])
    
    ## Get times
    ntj = nct$dim$time$len

    ## Load data
    mask = seq(t1, t1 + ntj - 1, 1)
    
    buffer = ncvar_get(nct, start = c(1,start,1), count = c(nx,count,ntj))
    temp0[,,mask] = buffer
    
    buffer = ncvar_get(ncp, start = c(1,start,1), count = c(nx,count,ntj))
    mask1 = which(buffer < 0)
    buffer[mask1] = 0
    mask1 = which(buffer < opts$threshold)
    buffer[mask1] = NA
    precip0[,,mask] = buffer
    
    if (exists("mask", opts)) {
      buffer = ncvar_get(ncm, start = c(1,start,1), count = c(nx,count,ntj))
      # mask0[,,mask] = !is.na(buffer) & buffer > 0
      mask0[,,mask] = buffer
    }    

    ## Close connections
    nc_close(nct)
    nc_close(ncp)
    if (exists("mask", opts))
      nc_close(ncm)
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

  ## Quantile regression
  print(paste("Quantile regression on chunk",i,"of",n.chunks))
  for (k in 1:nx) {
    for (l in 1:count) {
      temp1   = temp0  [k,l,]
      precip1 = precip0[k,l,]
      mask    = !is.na(precip1)
      if (exists("mask", opts)) {
        mask1 = mask0  [k,l,]
        mask  = mask & mask1
      }
      temp1   = temp1  [mask]
      precip1 = precip1[mask]
      qrm = rq(log(precip1) ~ temp1, tau = opts$quantile)
      qrm.summary = summary(qrm, se = "rank", alpha = 1 - opts$level)
      coefs = 100*(exp(qrm.summary$coefficients[2,]) - 1)
      beta [k,l] = coefs[1]
      lower[k,l] = coefs[2]
      upper[k,l] = coefs[3]
    } ## l
  } ## k
  rm(precip0,temp0,temp1,precip1,mask)
  if (exists("mask", opts))
    rm(mask0,mask1)
  gc()
  
  ## Transform data
  if (fliplon) {
    beta  = lonflip(beta , lon0)$x
    lower = lonflip(lower, lon0)$x
    upper = lonflip(upper, lon0)$x
  }
  if (fliplat & count > 1) {
    beta  = invertlat(beta )
    lower = invertlat(lower)
    upper = invertlat(upper)
  }
  
  ## Write data
  print(paste("Writing chunk",i,"of",n.chunks))
  if (fliplat) {
    ncvar_put(nco, "beta" , beta , 
              start = c(1,ny - start - count + 2), count = c(nx,count))
    ncvar_put(nco, "lower", lower, 
              start = c(1,ny - start - count + 2), count = c(nx,count))
    ncvar_put(nco, "upper", upper, 
              start = c(1,ny - start - count + 2), count = c(nx,count))
  } else {
    ncvar_put(nco, "beta" , beta , 
              start = c(1,start), count = c(nx,count))
    ncvar_put(nco, "lower", lower, 
              start = c(1,start), count = c(nx,count))
    ncvar_put(nco, "upper", upper, 
              start = c(1,start), count = c(nx,count))
  }

} ## i

## Close outout file
nc_close(nco)
