#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(ncdf4)
library(geosphere)

## Load source
source("src/datetime.R")
source("src/distance.R")
source("src/image.R")
source("src/invertlat.R")
source("src/lonflip.R")


filename = "erai_10fg_200001.nc"
varid = "fg10"

threshold = 250


## Open file
nci = nc_open(filename)
lon0 = as.numeric(nci$dim$lon$vals)
lat0 = as.numeric(nci$dim$lat$vals)
nlon = length(lon0)
nlat = length(lat0)
fliplon = lon0[1] < 0
fliplat = lat0[1] > lat0[2]
if (fliplon) {
  lon = fliplon(matrix(0, nlon, nlat), lon)
} else {
  lon = lon0
}
if (fliplat) {
  lat = rev(lat0)
} else {
  lat = lat0
}

times = ncvar_get(nci, "time")
time_units = ncatt_get(nci, "time", "units")$value
calendar = ncatt_get(nci, "time", "calendar")$value
facets = strsplit(time_units, " ")[[1]]
units = strsplit(time_units, " ")[[1]][1]
if (length(facets) == 3) {
  t0 = strptime(facets[3], format = "%Y-%m-%d", tz = "UTC")
} else if (length(facets) == 4) {
  t0 = strptime(paste(facets[3:4], collapse = " "),
                format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
} else if (length(facets) == 5) {
  t0 = strptime(paste(facets[3:5], collapse = " "),
                format = "%Y-%m-%d %H:%M:%OS %z", tz = "UTC")
}
t0 = format(t0, "%Y%m%dT%H%M%S")
dates = extract.dates(t0, times, calendar, units)
tt = length(times)

t = 1

input = ncvar_get(nci, varid, c(1,1,t), c(-1,-1,1))

if (fliplon) {
  buffer = lonflip(input, lon0)$x
  input = buffer$x
} else {
  lon = lon0
}
if (fliplat) {
  input = invertlat(input)
}


dd = distance(c(0,0), c(0,1))
dy = (lat[2] - lat[1])*dd

dny = floor(threshold/dy)


image(lon,lat,input, asp = 1)

# buffer = rbind(input,input)
buffer = rbind(input,input,input)

## Expand grid to take care of edge cases
xx = matrix(0, nrow(buffer) + 2, ncol(buffer) + 2)
xx[2:(nrow(buffer)+1),2:(ncol(buffer)+1)] = buffer

nx = nrow(xx)
ny = ncol(xx)

## Identify events
nevents = 0
events = list()
points = which(xx > 0)
npoints = length(points)
while(npoints > 0) {
  
  ## Increment counter
  nevents = nevents + 1
  
  ## Starting point
  x0 = points[npoints] %% nx
  y0 = points[npoints] %/% nx + 1
  
  ## Find all points connected to the starting point
  current = find.connected.component(xx, x0, y0)
  event = list()
  event$points = which(current > 0)
  event$x = event$points %% nx - 1
  event$y = event$points %/% nx #+ 1 - 1
  
  ## Store event
  events[[nevents]] = event
  
  ## Remove found points from master list
  mask = which(points %in% event$points)
  points = points[-mask]
  npoints = length(points)
  
}

## Remove events not partially in central region
i = 1
while(i <= length(events)) {
  # if (all(events[[i]][,1] > nlon)) {
  if (all(events[[i]]$x > 2*nlon) | any(events[[i]]$x <= nlon)) {
    events = events[-i]
  } else {
    # events[[i]]$x = events[[i]]$x - nlon
    events[[i]]$x = (events[[i]]$x - 1) %% nlon + 1
    events[[i]]$points = (events[[i]]$y - 1)*nlon + events[[i]]$x
    i = i + 1
  }
}
nevents = length(events)

output = matrix(NA, nlon, nlat)
for (i in 1:nevents)
  output[events[[i]]$points] = i

image(lon, lat, output, add = TRUE)
