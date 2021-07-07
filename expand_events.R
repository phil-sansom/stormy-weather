#!/usr/bin/env -S Rscript --vanilla

## Load libraries
library(ncdf4)

## Load source
source("src/datetime.R")
source("src/distance.R")
source("src/distances.R")
source("src/image.R")
source("src/invertlat.R")
source("src/lonflip.R")


filename = "../cyclone_id/tests/output/cyclones_20000102T00.nc"
varid = "idclust"

threshold = 250


## Open file
nci = nc_open(filename)
lon = as.numeric(nci$dim$lon$vals)
lat = as.numeric(nci$dim$lat$vals)
nlon = length(lon)
nlat = length(lat)

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

dd = distance(c(0,0), c(0,1))
dy = (lat[2] - lat[1])*dd

dny = floor(threshold/dy)

#input[input == 0] = NA
image(lon,lat,input, asp = 1)

ids = sort(unique(as.numeric(input)))
ids = ids[ids > 0 & !is.na(ids)]
nids = length(ids)

output = array(NA, c(nlon,nlat,nids))

for (i in 1:nids) {

  points  = which(input == ids[i])
  indices = arrayInd(points, c(nlon,nlat))

  for (j in 1:length(points)) {
    x = indices[j,1]
    y = indices[j,2]

    y1 = min(y + dny,nlat)
    y0 = max(y - dny, 1)

    dx  = (lon[2] - lon[1])*min(cos(pi*lat[y1]/180),cos(pi*lat[y0]/180))*dd
    dnx = floor(threshold/dx)

    x1 = min(x + dnx,nlon)
    x0 = max(x - dnx,1)

    grid = as.matrix(expand.grid(x0:x1,y0:y1))
    dists = distances(c(lon[x],lat[y]), cbind(lon[grid[,1]],lat[grid[,2]]))

    stencil = grid[dists <= threshold,]
    # stencil[,1] = stencil[,1] - x
    # stencil[,2] = stencil[,2] - y

    output[cbind(stencil,i)] = ids[i]

  } ## j

} ## i

check = apply(output, c(1,2), function(x) sum(!is.na(x)))

image(lon, lat, check)
contour(lon, lat, input, add = TRUE, levels = 1:max(ids))



# ## Bounding box in grid coordinates  
# minx = min(indices[,1])
# maxx = max(indices[,1])
# miny = min(indices[,2])
# maxy = max(indices[,2])
# 
# rect(lon[minx], lat[miny], lon[maxx], lat[maxy], lwd = 2)
# 
# ## Expand bounding box
# maxy1 = min(maxy + dny,nlat)
# miny1 = max(miny - dny, 1)
# 
# dx  = (lon[2] - lon[1])*min(cos(pi*lat[maxy1]/180),cos(pi*lat[miny1]/180))*dd
# dnx = floor(threshold/dx)
# 
# maxx1 = min(maxx + dnx,nlon)
# minx1 = max(minx - dnx,1)
# 
# rect(lon[minx1], lat[miny1], lon[maxx1], lat[maxy1], lwd = 2, border = "red")
