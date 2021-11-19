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
  make_option(c("--method","-m"), action = "store", type = "character",
              help = "Method of inference (Wald or rank) [default: Wald]",
              default = "Wald"),
  make_option(c("--nid","-n"), action = "store_false", type = "logical",
              help = "Not independent and identically distributed",
              dest = "iid", default = TRUE),
  make_option(c("--level","-l"), action = "store", type = "double",
              help = "Nominal confidence level required [default: 0.95]",
              default = 0.95),
  make_option("--compression", action = "store", type = "integer",
              help = "Compression level to use (0-9) [default: 5]",
              default = 5),
  make_option("--memory", action = "store", type = "integer",
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
if (!opts$method %in% c("Wald","rank"))
  stop("Method should be either Wald or rank")
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

## Precipitation
nc = nc_open(plist[1])
precip.name  = nc$var[[1]]$name
precip.units = nc$var[[1]]$units
precip.longname = nc$var[[1]]$longname
nc_close(nc)

## Mask levels
if (exists("mask", opts)) {
  nc = nc_open(mlist[1])
  mask.name = nc$var[[1]]$name
  levels = ncatt_get(nc, mask.name, "flag_values")
  labels = ncatt_get(nc, mask.name, "flag_meanings")
  if (levels$hasatt & labels$hasatt) {
    flags = TRUE
    levels = levels$value
    labels = strsplit(labels$value, " ")[[1]]
    nlevels = length(levels)
  } else {
    flags = FALSE
  }
}

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

## Parameter names
if (flags) {
  npar = 2*nlevels
  par.names = c("Intercept",labels[2:nlevels],
                temp.name,paste0(temp.name,":",labels[2:nlevels]))
} else {
  npar = 2
  par.names = c("Intercept",temp.name)
}
nchar = max(nchar(par.names))

## Define dimensions
lon.dim  = ncdim_def("longitude", "degrees_east" , lon, longname = "Longitude")
lat.dim  = ncdim_def("latitude" , "degrees_north", lat, longname = "Latitude")
par.dim  = ncdim_def("par"  , "", 1:npar , create_dimvar = FALSE)
char.dim = ncdim_def("nchar", "", 1:nchar, create_dimvar = FALSE)

## Define variables
par.var  = ncvar_def("parameter", "", list(char.dim,par.dim),
                     longname = "Parameter", prec = "char")
coef.var = ncvar_def("coefficients", "", list(lon.dim,lat.dim,par.dim),
                     missval = 9.9692099683868690e+36,
                     longname = "Coefficients", 
                     prec = "double",
                     compression = opts$compression)
vars = list(par.var,coef.var)

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

if (flags) {
  
  ## Define dimensions
  level.dim = ncdim_def("lev" , "", 1:nlevels, create_dimvar = FALSE)
  test.dim  = ncdim_def("test", "", 1:3      , create_dimvar = FALSE)

  ## Define variables
  level.var = ncvar_def("level", "", list(char.dim,level.dim),
                        longname = "Level", prec = "char")
  test.var  = ncvar_def("ftest", "", list(char.dim,test.dim),
                        longname = "F-test", prec = "char")
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
  vars = c(vars, list(level.var, test.var, count.var, pval.var))
  
} ## flags

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
ncatt_put(nco, 0, "method"   , opts$method   , prec = "text"   )
ncatt_put(nco, 0, "iid"      , opts$iid      , prec = "integer")
if (opts$method == "rank")
  ncatt_put(nco, 0, "level", opts$level, prec = "double")

## Write auxiliary coordinate variable
ncvar_put(nco, "parameter", par.names)
ncatt_put(nco, "coefficients", "coordinates", "parameter")
if (opts$method == "Wald") {
  ncatt_put(nco, "covariance", "coordinates", "parameter parameter")
} else if (opts$method == "rank") {
  ncatt_put(nco, "lower", "coordinates", "parameter")
  ncatt_put(nco, "upper", "coordinates", "parameter")
}
if (flags) {
  ncvar_put(nco, "level", labels)
  ncvar_put(nco, "ftest", c(mask.name,temp.name,"interaction"))
  ncatt_put(nco, "counts", "coordinates", "level")
  ncatt_put(nco, "pvalue", "coordinates", "ftest")
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
  if (exists("mask", opts))
    mask0 = array(NA, c(nx,count,nt))
  if (flags) {
    counts  = array(NA, c(nx,count,nlevels))
    pvalues = array(NA, c(nx,count,3))
  }
  coef = array(NA, c(nx,count,npar))
  if (opts$method == "Wald") {
    cov   = array(NA, c(nx,count,npar,npar))
    df    = array(NA, c(nx,count))
  } else if (opts$method == "rank") {
    lower = array(NA, c(nx,count,npar))
    upper = array(NA, c(nx,count,npar))
  }

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
      if (flags) {
        mask0[,,mask] = buffer
      } else {
        mask0[,,mask] = !is.na(buffer) & buffer > 0
      }
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
    
  } ## smoothing

  ## Quantile regression
  print(paste("Quantile regression on chunk",i,"of",n.chunks))
  for (k in 1:nx) {
    for (l in 1:count) {
      
      ## Extract data
      temp1   = temp0  [k,l,]
      precip1 = precip0[k,l,]
      mask    = !is.na(precip1)
      if (exists("mask", opts)) {
        mask1 = mask0  [k,l,]
        if (flags) {
          mask1 = mask1[mask]
          mask1 = factor(mask1, levels, labels)
        } else {
          mask  = mask & mask1
        }
      }
      temp1   = temp1  [mask]
      precip1 = precip1[mask]

      ## Fit model
      if (flags) {

        counts[k,l,] = as.numeric(table(mask1))
        
        rqm = try(rq(log(precip1) ~ mask1 + temp1 + mask1:temp1, 
                     tau = opts$quantile), TRUE)
        
      } else {
        
        rqm = try(rq(log(precip1) ~ temp1, 
                     tau = opts$quantile, iid = opts$iid), TRUE)

      } ## flags
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
      
      ## Analysis of deviance
      if (flags) {
        
        ## Fit simpler models for comparison
        rq0 = try(rq(log(precip1) ~ 1,
                     tau = opts$quantile), TRUE)
        if (class(rq0) == "try-error")
          next
        rq1 = try(rq(log(precip1) ~ mask1,
                     tau = opts$quantile), TRUE)
        if (class(rq1) == "try-error")
          next
        rq2 = try(rq(log(precip1) ~ mask1 + temp1, 
                     tau = opts$quantile), TRUE)
        if (class(rq2) == "try-error")
          next
        
        ## Compute analysis of deviance
        aod = try(anova(rqm, rq2, rq1, rq0, test = opts$method, 
                    se = ifelse(opts$iid, "iid", "nid"), iid = opts$iid), TRUE)
        if (class(aod) == "try-error")
          next

        ## Store p-values
        pvalues[k,l,] = rev(aod$table[,4])
        
      } ## flags
      
      ## Extract coefficients
      coefficients = rqm$coefficients
      
      ## Check that all levels populated
      if (length(coefficients) < npar) {
        
        buffer = coefficients
        coefficients = rep(NA, npar)
        missing = which(! labels %in% rqm$xlevels$mask1)
        slice = rep(TRUE, npar)
        slice[c(missing,npar/2+missing)] = FALSE
        coefficients[slice] = buffer
        
      }
      
      ## Store coefficients
      coef[k,l,] = coefficients

      if (opts$method == "Wald") {
        
        ## Extract covariance
        covariance = rqm.summary$cov
        
        ## Check that all levels populated
        if (nrow(covariance) < npar) {
          
          buffer = covariance
          covariance = matrix(NA, npar, npar)
          missing = which(! labels %in% rqm$xlevels$mask1)
          slice = matrix(TRUE, npar, npar)
          slice[c(missing,npar/2+missing),] = FALSE
          slice[,c(missing,npar/2+missing)] = FALSE
          covariance[slice] = buffer

        }
        
        ## Store results
        cov[k,l,,] = covariance
        df [k,l  ] = rqm.summary$rdf
        
      } else if (opts$method == "rank") {
        
        ## Extract bounds
        bounds = rqm.summary$coefficients[,2:3]
        
        ## Check that all levels populated
        if (nrow(bounds) < npar) {

          buffer = bounds
          bounds = matrix(NA, npar, 2)
          missing = which(! labels %in% rqm$xlevels$mask1)
          slice = matrix(TRUE, npar, 2)
          slice[c(missing,npar/2+missing),] = FALSE
          bounds[slice] = buffer
          
        }
        
        ## Store results
        lower[k,l,] = bounds[,1]
        upper[k,l,] = bounds[,2]
      
      }
      
      
    } ## l
  } ## k
  rm(precip0,temp0,temp1,precip1,mask)
  if (exists("mask", opts))
    rm(mask0,mask1)
  gc()
  
  ## Transform data
  if (fliplon) {
    coef  = lonflip(coef , lon0)$x
    if (opts$method == "Wald") {
      cov   = lonflip(cov  , lon0)$x
      df    = lonflip(df   , lon0)$x
    } else if (opts$method == "rank") {
      lower = lonflip(lower, lon0)$x
      upper = lonflip(upper, lon0)$x
    }
    if (flags) {
      counts = lonflip(counts, lon0)$x
    }
  }
  if (fliplat & count > 1) {
    coef  = invertlat(coef )
    if (opts$method == "Wald") {
      cov   = invertlat(cov)
      df    = invertlat(df )
    } else if (opts$method == "rank") {
      lower = invertlat(lower)
      upper = invertlat(upper)
    }
    if (flags) {
      counts = invertlat(counts)
    }
  }
  
  ## Write data
  print(paste("Writing chunk",i,"of",n.chunks))
  if (fliplat) {
    ncvar_put(nco, "coefficients", coef, 
              start = c(1,ny - start - count + 2,1), 
              count = c(nx,count,npar))
    if (opts$method == "Wald") {
      ncvar_put(nco, "covariance", cov, 
                start = c(1,ny - start - count + 2,1,1), 
                count = c(nx,count,npar,npar))
      ncvar_put(nco, "df", df, 
                start = c(1,ny - start - count + 2),
                count = c(nx,count))
    } else if (opts$method == "rank") {
      ncvar_put(nco, "lower", lower, 
                start = c(1,ny - start - count + 2,1), 
                count = c(nx,count,npar))
      ncvar_put(nco, "upper", upper, 
                start = c(1,ny - start - count + 2,1),
                count = c(nx,count,npar))
    }
    if (flags) {
      ncvar_put(nco, "counts", counts, 
                start = c(1,ny - start - count + 2,1), 
                count = c(nx,count,nlevels))
    }
  } else {
    ncvar_put(nco, "coefficients", coef, 
              start = c(1,start,1), count = c(nx,count,npar))
    if (opts$method == "Wald") {
      ncvar_put(nco, "covariance", cov, 
                start = c(1,start,1,1), count = c(nx,count,npar,npar))
      ncvar_put(nco, "df", df, 
                start = c(1,start), count = c(nx,count))
    } else if (opts$method == "rank") {
      ncvar_put(nco, "lower", lower, 
                start = c(1,start,1), count = c(nx,count,npar))
      ncvar_put(nco, "upper", upper, 
                start = c(1,start,1), count = c(nx,count,npar))
    }
    if (flags) {
      ncvar_put(nco, "counts", counts, 
                start = c(1,start,1), 
                count = c(nx,count,nlevels))
    }
  }

} ## i

## Close outout file
nc_close(nco)
