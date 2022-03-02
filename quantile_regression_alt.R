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
  make_option(c("--threshold","-t"), action = "store", type = "double",
              help = "Threshold for wet days (mm) [default: 0.1]",
              default = 0.1),
  make_option(c("--method","-m"), action = "store", type = "character",
              help = "Method of inference (Wald or rank) [default: Wald]",
              default = "Wald"),
  make_option(c("--nid","-n"), action = "store_false", type = "logical",
              help = "Not independent and identically distributed",
              dest = "iid", default = TRUE),
  make_option(c("--level","-l"), action = "store", type = "double",
              help = "Nominal confidence level required [default: 0.95]",
              default = 0.95),
  make_option(c("--none"), action = "store_true", type = "logical",
              help = "Include the no storm type",
              default = "FALSE"),
  make_option("--compression", action = "store", type = "integer",
              help = "Compression level to use (0-9) [default: 5]",
              default = 5L),
  make_option("--memory", action = "store", type = "integer",
              help = "Maximum memory to use (in MB)")
)

## Argument parser
parser = OptionParser(usage = "Usage: %prog [OPTION]... CYCLONE_FILES FRONT_FILES THUNDER_FILES TEMP_FILES PRECIP_FILES OUTFILE",
                      option_list = option_list,
                      description = "Clausius-Clapeyron scaling calculations.\n\nOperands:\n\tCYCLONE_FILES\n\t\tList of cyclone footprint files to read\n\tFRONT_FILES\n\t\tList of front footprint files to read\n\tTHUNDER_FILES\n\t\tList of thunder footprint files to read\n\tTEMP_FILES\n\t\tList of temperature files to read\n\tPRECIP_FILES\n\t\tList of precipitation files to read\n\tOUTFILE\n\t\tFile to write output to")

## Parser arguments
args = commandArgs(trailingOnly = TRUE)
argv = parse_args(parser, args, positional_arguments = 6L)
opts = argv$options
args = argv$args
if (!opts$method %in% c("Wald","rank"))
  stop("Method should be either Wald or rank")
if (opts$compression == 0L)
  opts$compression = NA

## Read file lists
clist = scan(args[1L], character(), -1L, quiet = TRUE)
flist = scan(args[2L], character(), -1L, quiet = TRUE)
tlist = scan(args[3L], character(), -1L, quiet = TRUE)
dlist = scan(args[4L], character(), -1L, quiet = TRUE)
plist = scan(args[5L], character(), -1L, quiet = TRUE)
if (length(flist) != length(clist))
  stop("Length of front file list differs from length of cyclone file list")
if (length(tlist) != length(flist))
  stop("Length of thunder file list differs from length of front file list")
if (length(dlist) != length(tlist))
  stop("Length of temperature file list differs from length of thunder file list")
if (length(plist) != length(dlist))
  stop("Length of precipitation file list differs from length of temperature file list")
n.files = length(clist)

## Check dimensions
print("Checking dimensions...")
time = numeric()
for (i in 1L:n.files) {
  
  ## Cyclone file
  ncc = nc_open(clist[i])
  lonc = ncc$dim$longitude$vals
  latc = ncc$dim$latitude$vals
  timec = ncc$dim$time$vals
  nc_close(ncc)
  
  ## Front file
  ncf = nc_open(flist[i])
  lonf = ncf$dim$longitude$vals
  latf = ncf$dim$latitude$vals
  timef = ncf$dim$time$vals
  nc_close(ncf)
  
  ## Thunder file
  nct = nc_open(tlist[i])
  lont = nct$dim$longitude$vals
  latt = nct$dim$latitude$vals
  timet = nct$dim$time$vals
  nc_close(nct)
  
  ## Temperature file
  ncd = nc_open(dlist[i])
  lond = ncd$dim$longitude$vals
  latd = ncd$dim$latitude$vals
  timed = ncd$dim$time$vals
  nc_close(ncd)
  
  ## Precipitation file
  ncp = nc_open(plist[i])
  lonp = ncp$dim$longitude$vals
  latp = ncp$dim$latitude$vals
  timep = ncp$dim$time$vals
  nc_close(ncc)
  
  ## Check dimensions
  if (any(lonf != lonc))
    stop(paste("Longitude dimensions don't match between front and cyclone files", i))
  if (any(lont != lonf))
    stop(paste("Longitude dimensions don't match between thunder and front files", i))
  if (any(lond != lont))
    stop(paste("Longitude dimensions don't match between temperature and thunder files", i))
  if (any(lonp != lond))
    stop(paste("Longitude dimensions don't match between precipitation and temperature files", i))
  if (any(latf != latc))
    stop(paste("Latitude dimensions don't match between front and cyclone files", i))
  if (any(latt != latf))
    stop(paste("Latitude dimensions don't match between thunder and front files", i))
  if (any(latd != latt))
    stop(paste("Latitude dimensions don't match between temperature and thunder files", i))
  if (any(latp != latd))
    stop(paste("Latitude dimensions don't match between precipitation and temperature files", i))
  if (any(timef != timec))
    stop(paste("Time dimensions don't match between front and cyclone files", i))
  if (any(timet != timef))
    stop(paste("Time dimensions don't match between thunder and front files", i))
  if (any(timed != timet))
    stop(paste("Time dimensions don't match between temperature and thunder files", i))
  if (any(timep != timed))
    stop(paste("Time dimensions don't match between precipitation and temperature files", i))
  
  if (i > 1) {
    if (any(lonc != lonc0))
      stop(paste("Longitude dimensions don't match between cyclone files", i-1L, "and", i))
    if (any(latc != latc0))
      stop(paste("Latitude dimensions don't match between cyclone files",  i-1L, "and", i))
    if (any(lonf != lonf0))
      stop(paste("Longitude dimensions don't match between front files", i-1L, "and", i))
    if (any(latf != latf0))
      stop(paste("Latitude dimensions don't match between front files",  i-1L, "and", i))
    if (any(lont != lont0))
      stop(paste("Longitude dimensions don't match between thunder files", i-1L, "and", i))
    if (any(latt != latt0))
      stop(paste("Latitude dimensions don't match between thunder files",  i-1L, "and", i))
    if (any(lond != lond0))
      stop(paste("Longitude dimensions don't match between temperature files", i-1L, "and", i))
    if (any(latd != latd0))
      stop(paste("Latitude dimensions don't match between temperature files",  i-1L, "and", i))
    if (any(lonp != lonp0))
      stop(paste("Longitude dimensions don't match between precipitation files", i-1L, "and", i))
    if (any(latp != latp0))
      stop(paste("Latitude dimensions don't match between precipitation files",  i-1L, "and", i))
  }
  lonc0 = lonc; latc0 = latc; timec0 = timec
  lonf0 = lonf; latf0 = latf; timef0 = timef
  lont0 = lont; latt0 = latt; timet0 = timet
  lond0 = lond; latd0 = latd; timed0 = timed
  lonp0 = lonp; latp0 = latp; timep0 = timep
  
  time = c(time,timec)
  
} ## i
nt = length(time)
dt = time[2L] - time[1L]

## Read dimensions
nc = nc_open(dlist[1])
lon0 = nc$dim$longitude$vals
lat0 = nc$dim$latitude$vals
calendar   = nc$dim$time$calendar
time.units = nc$dim$time$units
temp.name  = nc$var[[1L]]$name
temp.units = nc$var[[1L]]$units
temp.longname = nc$var[[1L]]$longname
nc_close(nc)

## Precipitation
nc = nc_open(plist[1L])
precip.name  = nc$var[[1L]]$name
precip.units = nc$var[[1L]]$units
precip.longname = nc$var[[1L]]$longname
nc_close(nc)

## Labels
if (opts$none) {
levels = c("none","cyclone","front","thunder",
           "cyclone:front","cyclone:thunder","front:thunder",
           "cyclone:front:thunder")
labels = c("(Intercept)","cyclone","front","thunder","temperature",
           "cyclone:front","cyclone:thunder","front:thunder",
           "cyclone:temperature","front:temperature","thunder:temperature",
           "cyclone:front:thunder","cyclone:front:temperature",
           "cyclone:thunder:temperature","front:thunder:temperature",
           "cyclone:front:thunder:temperature")
} else {
  levels = c("cyclone","front","thunder",
             "cyclone:front","cyclone:thunder","front:thunder",
             "cyclone:front:thunder")
  labels = c("(Intercept)","front","thunder","temperature",
             "cyclone:front","cyclone:thunder","front:thunder",
             "front:temperature","thunder:temperature",
             "cyclone:front:thunder","cyclone:front:temperature",
             "cyclone:thunder:temperature","front:thunder:temperature",
             "cyclone:front:thunder:temperature")
}
tests   = c("storm_type","temperature","storm_type:temperature")
nlevels = length(levels)
nlabels = length(labels)
ntests  = length(tests)

if (precip.units %in% c("mm","cm","m")) {
  intercept.units = "ln(mm)"
  slope.units     = paste0("ln(mm)/",temp.units)
} else {
  intercept.units = paste0("ln(",precip.units,")")
  slope.units     = paste0("ln(",precip.units,")/",temp.units)
}

nx = length(lon0)
ny = length(lat0)

## Determine if we need to transform dimensions
fliplon = any(lon0 >= 180)
fliplat = lat0[1L] > lat0[2L]

## Transform dimensions
if (fliplon) {
  lon = lonflip(matrix(0, length(lon0), length(lat0)), lon0)$lon 
} else {
  lon = lon0
} ## fliplon
lat = if (fliplat) rev(lat0) else lat0

## Split into chunks
if (exists("memory", opts)) {
  row.size   = nx*nt*8/1024/1024
  chunk.size = floor(opts$memory/row.size/5)
  n.chunks   = ceiling(ny/chunk.size)
} else {
  chunk.size = ny
  n.chunks = 1L
} ## memory
chunks = data.frame(
  start = seq(0L, n.chunks - 1L, 1L)*chunk.size + 1L,
  count = rep(chunk.size, n.chunks)
)
chunks$count[n.chunks] = ny - (n.chunks - 1L)*chunk.size
if (fliplat) {
  chunks$start = rev(chunks$start)
  chunks$count = rev(chunks$count)
} ## fliplat

## Parameter names
npar = nlabels
par.names = labels
nchar = max(nchar(par.names),nchar(tests))

## Define dimensions
lon.dim   = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim   = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude")
par.dim   = ncdim_def("par"  , "", 1L:npar   , create_dimvar = FALSE)
level.dim = ncdim_def("lev"  , "", 1L:nlevels, create_dimvar = FALSE)
test.dim  = ncdim_def("test" , "", 1L:ntests , create_dimvar = FALSE)
char.dim  = ncdim_def("nchar", "", 1L:nchar  , create_dimvar = FALSE)

## Define variables
level.var = ncvar_def("level", "", list(char.dim,level.dim),
                      longname = "Level", prec = "char")
par.var   = ncvar_def("parameter", "", list(char.dim,par.dim),
                      longname = "Parameter", prec = "char")
test.var  = ncvar_def("ftest", "", list(char.dim,test.dim),
                      longname = "F-test", prec = "char")
vars = list(level.var,par.var,test.var)
count.var = ncvar_def("counts", "", list(lon.dim,lat.dim,level.dim),
                      missval = -2147483647L,
                      longname = "Counts", 
                      prec = "integer",
                      compression = opts$compression)
pval.var = ncvar_def("pvalue", "", list(lon.dim,lat.dim,test.dim),
                     missval = 9.9692099683868690e+36,
                     longname = "p-value", 
                     prec = "double",
                     compression = opts$compression)
coef.var = ncvar_def("coefficients", "", list(lon.dim,lat.dim,par.dim),
                     missval = 9.9692099683868690e+36,
                     longname = "Coefficients", 
                     prec = "double",
                     compression = opts$compression)
vars = c(vars,list(count.var,pval.var,coef.var))

if (opts$method == "Wald") {
  
  ## Define variables
  cov.var = ncvar_def("covariance", "", list(lon.dim,lat.dim,par.dim,par.dim),
                      missval = 9.9692099683868690e+36,
                      longname = "Covariance matrix", 
                      prec = "double",
                      compression = opts$compression)
  df.var  = ncvar_def("df", "", list(lon.dim,lat.dim),
                      missval = -2147483647L,
                      longname = "Residual degrees of freedom", 
                      prec = "integer",
                      compression = opts$compression)
  vars = c(vars, list(cov.var, df.var))
  
} else if (opts$method == "rank") {
  
  ## Define variables
  lwr.var = ncvar_def("lower", "", list(lon.dim,lat.dim,par.dim),
                      missval = 9.9692099683868690e+36,
                      longname = "Lower bound", 
                      prec = "double",
                      compression = opts$compression)
  upr.var = ncvar_def("upper", "", list(lon.dim,lat.dim,par.dim),
                      missval = 9.9692099683868690e+36,
                      longname = "Upper bound", 
                      prec = "double",
                      compression = opts$compression)
  vars = c(vars, list(lwr.var, upr.var))
  
}

## Create netCDF file
nco = nc_create(args[6L], vars)

## Write description
ncatt_put(nco, 0L, "Conventions", "CF-1.9", prec = "text")

## Write standard names
ncatt_put(nco, "longitude", "standard_name", "longitude", prec = "text")
ncatt_put(nco, "latitude" , "standard_name", "latitude" , prec = "text")

## Write axes
ncatt_put(nco, "longitude", "axis", "X", prec = "text")
ncatt_put(nco, "latitude" , "axis", "Y", prec = "text")

## Write parameters
ncatt_put(nco, 0L, "quantile" , opts$quantile , prec = "double" )
ncatt_put(nco, 0L, "threshold", opts$threshold, prec = "double" )
ncatt_put(nco, 0L, "none"     , opts$none     , prec = "integer")
ncatt_put(nco, 0L, "method"   , opts$method   , prec = "text"   )
ncatt_put(nco, 0L, "iid"      , opts$iid      , prec = "integer")
if (opts$method == "rank")
  ncatt_put(nco, 0L, "level", opts$level, prec = "double")

## Write auxiliary coordinate variable
ncvar_put(nco, "parameter", par.names)
ncatt_put(nco, "coefficients", "coordinates", "parameter")
if (opts$method == "Wald") {
  ncatt_put(nco, "covariance", "coordinates", "parameter parameter")
} else if (opts$method == "rank") {
  ncatt_put(nco, "lower", "coordinates", "parameter")
  ncatt_put(nco, "upper", "coordinates", "parameter")
} ## method
ncvar_put(nco, "level", levels)
ncvar_put(nco, "ftest", tests )
ncatt_put(nco, "counts", "coordinates", "level")
ncatt_put(nco, "pvalue", "coordinates", "ftest")

## Transform threshold
if (precip.units == "m") {
  opts$threshold = opts$threshold/1000
} else if (precip.units == "cm") {
  opts$threshold = opts$threshold/10
} else if (! precip.units == "mm") {
  warning("Precip units not recognised (mm,cm,m), specify threshold in native units")
} ## precip.units

## Loop over chunks
for (i in 1L:n.chunks) {
  
  ## Initialize storage
  start = chunks$start[i]
  count = chunks$count[i]
  cyclone0 = array(NA, c(nx,count,nt))
  front0   = array(NA, c(nx,count,nt))
  thunder0 = array(NA, c(nx,count,nt))
  temp0    = array(NA, c(nx,count,nt))
  precip0  = array(NA, c(nx,count,nt))
  counts  = array(NA, c(nx,count,nlevels))
  pvalues = array(NA, c(nx,count,ntests))
  coef = array(NA, c(nx,count,npar))
  if (opts$method == "Wald") {
    cov   = array(NA, c(nx,count,npar,npar))
    df    = array(NA, c(nx,count))
  } else if (opts$method == "rank") {
    lower = array(NA, c(nx,count,npar))
    upper = array(NA, c(nx,count,npar))
  }
  
  ## Initialize time counter
  t1 = 1L
  
  ## Loop over files
  for (j in 1L:n.files) {
    
    ## Print status
    print(paste0("Reading chunk ",i," of ",n.chunks,
                 ", file ",j," of ", n.files))
    
    ## Open connections
    ncc = nc_open(clist[j])
    ncf = nc_open(flist[j])
    nct = nc_open(tlist[j])
    ncd = nc_open(dlist[j])
    ncp = nc_open(plist[j])

    ## Get times
    ntj = ncd$dim$time$len
    
    ## Time storage mask
    mask = seq(t1, t1 + ntj - 1L, 1L)
    
    ## Read cyclones
    buffer = ncvar_get(ncc, start = c(1L,start,1L), count = c(nx,count,ntj))
    cyclone0[,,mask] = buffer
    
    ## Read fronts
    buffer = ncvar_get(ncf, start = c(1L,start,1L), count = c(nx,count,ntj))
    front0[,,mask] = buffer
    
    ## Read thunder
    buffer = ncvar_get(nct, start = c(1L,start,1L), count = c(nx,count,ntj))
    thunder0[,,mask] = buffer
    
    ## Read temperature
    buffer = ncvar_get(ncd, start = c(1L,start,1L), count = c(nx,count,ntj))
    temp0[,,mask] = buffer
    
    ## Read precipitation
    buffer = ncvar_get(ncp, start = c(1L,start,1L), count = c(nx,count,ntj))
    mask1 = which(buffer < 0)
    buffer[mask1] = 0
    mask1 = which(buffer < opts$threshold)
    buffer[mask1] = NA
    precip0[,,mask] = buffer
    
    ## Close connections
    nc_close(ncc)
    nc_close(ncf)
    nc_close(nct)
    nc_close(ncd)
    nc_close(ncp)
    
    ## Garbage collection
    rm(buffer,mask,mask1)
    gc()
    
    ## Increment time counter
    t1 = t1 + ntj
    
  } ## j
  gc()
  
  
  #########################
  ## Quantile regression ##
  #########################
  
  print(paste("Quantile regression on chunk",i,"of",n.chunks))
  
  ## Loop over longitudes
  for (k in 1L:nx) {

    ## Loop over latitudes
    for (l in 1L:count) {
      
      ## Extract data
      precipitation = precip0 [k,l,]
      cyclone       = cyclone0[k,l,]
      front         = front0  [k,l,]
      thunder       = thunder0[k,l,]
      temperature   = temp0   [k,l,]
      
      if (opts$none) {
        mask = !is.na(precipitation)
      } else {
        mask = !is.na(precipitation) & (cyclone | front | thunder)
      }
        
      precipitation = precipitation[mask]
      cyclone       = cyclone      [mask]
      front         = front        [mask]
      thunder       = thunder      [mask]
      temperature   = temperature  [mask]

      ## Counts
      nC   = sum(cyclone)
      nF   = sum(front)
      nT   = sum(thunder)
      nCF  = sum(cyclone & front)
      nCT  = sum(cyclone & thunder)
      nFT  = sum(front   & thunder)
      nCFT = sum(cyclone & front & thunder)
      if (opts$none) {
        nN = sum(!(cyclone | front | thunder))
        counts[k,l,] = c(nN,nC,nF,nT,nCF,nCT,nFT,nCFT)
      } else {
        counts[k,l,] = c(nC,nF,nT,nCF,nCT,nFT,nCFT)
      }
      
      ## Fit model
      if (opts$none) {
        rqm = try(rq(log(precipitation) ~ cyclone*front*thunder*temperature, 
                     tau = opts$quantile, iid = opts$iid), TRUE)
      } else {
        rqm = try(rq(log(precipitation) ~ front + thunder + 
                       cyclone:front + cyclone:thunder + front:thunder + 
                       cyclone:front:thunder + temperature + 
                       front:temperature + thunder:temperature + 
                       cyclone:front:temperature + cyclone:thunder:temperature +
                       front:thunder:temperature + 
                       cyclone:front:thunder:temperature, 
                     tau = opts$quantile, iid = opts$iid), TRUE)
      }
      if (class(rqm) == "try-error")
        next
      
      ## Model summary
      if (opts$method == "Wald") {
        
        rqm.summary = try(summary(rqm, se = ifelse(opts$iid, "iid", "nid"), 
                                  covariance = TRUE), TRUE)
        
      } else if (opts$method == "rank") {
        
        rqm.summary = try(summary(rqm, se = "rank", alpha = 1 - opts$level, 
                                  iid = opts$iid), TRUE)
        
      } ## method
      if (class(rqm.summary) == "try-error")
        next
      
      ## Fit null model
      rq0 = try(rq(log(precipitation) ~ 1L, 
                   tau = opts$quantile, iid = opts$iid), TRUE)
      if (class(rq0) == "try-error")
        next
      
      ## Fit intercepts
      if (opts$none) {
        rq1 = try(rq(log(precipitation) ~ cyclone*front*thunder, 
                     tau = opts$quantile, iid = opts$iid), TRUE)
      } else {
        rq1 = try(rq(log(precipitation) ~ front + thunder + 
                       cyclone:front + cyclone:thunder + front:thunder + 
                       cyclone:front:thunder,
                     tau = opts$quantile, iid = opts$iid), TRUE)
      }
      if (class(rq0) == "try-error")
        next

      ## Fit slope
      if (opts$none) {
        rq2 = try(rq(log(precipitation) ~ cyclone*front*thunder + temperature, 
                     tau = opts$quantile, iid = opts$iid), TRUE)
      } else {
        rq2 = try(rq(log(precipitation) ~ front + thunder + 
                       cyclone:front + cyclone:thunder + front:thunder + 
                       cyclone:front:thunder + temperature,
                     tau = opts$quantile, iid = opts$iid), TRUE)
      }
      if (class(rq0) == "try-error")
        next
      
      ## Compute analysis of deviance
      aod = try(anova(rqm, rq2, rq1, rq0, test = opts$method, 
                      se = ifelse(opts$iid, "iid", "nid"), iid = opts$iid), 
                TRUE)
      if (class(aod) == "try-error")
        next
      
      ## Store p-values
      pvalues[k,l,] = rev(aod$table[,4L])
      
      ## Extract coefficients
      coefficients = rqm$coefficients
      
      ## Check that all levels populated
      if (length(coefficients) < npar) {
        
        present = which(par.names %in% names(coefficients))
        buffer = coefficients
        coefficients = rep(NA, npar)
        coefficients[present] = buffer
        
      }
      
      ## Store coefficients
      coef[k,l,] = coefficients
      
      if (opts$method == "Wald") {
        
        ## Extract covariance
        covariance = rqm.summary$cov
        
        ## Check that all levels populated
        if (nrow(covariance) < npar) {
          
          present = which(par.names %in% names(coefficients))
          buffer = covariance
          covariance = matrix(NA, npar, npar)
          covariance[present,present] = buffer
          
        }
        
        ## Store results
        cov[k,l,,] = covariance
        df [k,l  ] = rqm.summary$rdf
        
      } else if (opts$method == "rank") {
        
        ## Extract bounds
        bounds = rqm.summary$coefficients[,2L:3L]
        
        ## Check that all levels populated
        if (nrow(bounds) < npar) {
          
          present = which(par.names %in% names(coefficients))
          buffer = bounds
          bounds = matrix(NA, npar, 2L)
          bounds[present,] = buffer
          
        }
        
        ## Store results
        lower[k,l,] = bounds[,1L]
        upper[k,l,] = bounds[,2L]
        
      } ## method
      
    } ## l
  } ## k
  rm(precip0,temp0,cyclone0,front0,thunder0,
     temperature,precipitation,cyclone,front,thunder,mask)
  gc()
  
  ## Transform data
  if (fliplon) {
    for (j in 1L:nlevels)
      counts[,,j] = lonflip(counts[,,j], lon0)$x
    for (j in 1L:ntests)
      pvalues[,,j] = lonflip(pvalues[,,j], lon0)$x
    for (j in 1L:npar)
      coef[,,j] = lonflip(coef[,,j], lon0)$x
    if (opts$method == "Wald") {
      for (j in 1L:npar)
        for (k in 1L:npar)
          cov[,,j,k] = lonflip(cov[,,j,k], lon0)$x
      df = lonflip(df, lon0)$x
    } else if (opts$method == "rank") {
      for (j in 1L:npar) {
        lower[,,j] = lonflip(lower[,,j], lon0)$x
        upper[,,j] = lonflip(upper[,,j], lon0)$x
      } ## npar
    } ## method
  } ## fliplon
  if (fliplat & count > 1L) {
    for (j in 1L:nlevels)
      counts[,,j] = invertlat(counts[,,j])
    for (j in 1L:ntests)
      pvalues[,,j] = invertlat(pvalues[,,j])
    for (j in 1L:npar)
      coef[,,j] = invertlat(coef[,,j])
    if (opts$method == "Wald") {
      for (j in 1L:npar)
        for (k in 1L:npar)
          cov[,,j,k] = invertlat(cov[,,j,k])
      df = invertlat(df )
    } else if (opts$method == "rank") {
      for (j in 1L:npar) {
        lower[,,j] = invertlat(lower[,,j])
        upper[,,j] = invertlat(upper[,,j])
      } ## npar
    } ## method
  } ## fliplat
  
  ## Write data
  print(paste("Writing chunk",i,"of",n.chunks))
  if (fliplat) {
    ncvar_put(nco, "counts", counts, 
              start = c(1L,ny - start - count + 2L,1L), 
              count = c(nx,count,nlevels))
    ncvar_put(nco, "pvalue", pvalues, 
              start = c(1L,ny - start - count + 2L,1L), 
              count = c(nx,count,ntests))
    ncvar_put(nco, "coefficients", coef, 
              start = c(1L,ny - start - count + 2L,1L), 
              count = c(nx,count,npar))
    if (opts$method == "Wald") {
      ncvar_put(nco, "covariance", cov, 
                start = c(1L,ny - start - count + 2L,1L,1L), 
                count = c(nx,count,npar,npar))
      ncvar_put(nco, "df", df, 
                start = c(1L,ny - start - count + 2L),
                count = c(nx,count))
    } else if (opts$method == "rank") {
      ncvar_put(nco, "lower", lower, 
                start = c(1L,ny - start - count + 2L,1L), 
                count = c(nx,count,npar))
      ncvar_put(nco, "upper", upper, 
                start = c(1L,ny - start - count + 2L,1L),
                count = c(nx,count,npar))
    }
  } else {
    ncvar_put(nco, "counts", counts, 
              start = c(1L,start,1L), count = c(nx,count,nlevels))
    ncvar_put(nco, "pvalue", pvalues, 
              start = c(1L,start,1L), count = c(nx,count,ntests))
    ncvar_put(nco, "coefficients", coef, 
              start = c(1L,start,1L), count = c(nx,count,npar))
    if (opts$method == "Wald") {
      ncvar_put(nco, "covariance", cov, 
                start = c(1L,start,1L,1L), count = c(nx,count,npar,npar))
      ncvar_put(nco, "df", df, 
                start = c(1L,start), count = c(nx,count))
    } else if (opts$method == "rank") {
      ncvar_put(nco, "lower", lower, 
                start = c(1L,start,1L), count = c(nx,count,npar))
      ncvar_put(nco, "upper", upper, 
                start = c(1L,start,1L), count = c(nx,count,npar))
    }
  } ## fliplat
  
} ## i

## Close output file
nc_close(nco)
