find.connected.component = function (z, x0, y0) {
  
  # https://en.wikipedia.org/wiki/Connected-component_labeling
  
  target = z[x0,y0]
  mask = matrix(0, nrow(z), ncol(z))
  queue = list()
  queue[[1]] = c(x0,y0)
  mask[x0,y0] = 1
  while (length(queue) > 0) {
    current = queue[[1]]
    queue = queue[-1]
    x = current[1]
    y = current[2]
    zz = z[(x-1):(x+1),(y-1):(y+1)]
    new = which(zz == target)
    for (n in new) {
      xn = x + (n-1) %% 3  - 1
      yn = y + (n-1) %/% 3 - 1 
      if (mask[xn,yn] == 0) {
        queue[[length(queue) + 1]] = c(xn,yn)
        mask[xn,yn] = 1
      }
    } ## n
  }
  
  return(mask)
  
}

is.connected.minima = function(z, x0, y0) {
  
  component = find.connected.component(z, x0, y0)
  indices = which(component == 1, arr.ind = TRUE)
  values = numeric(0)
  for (i in 1:nrow(indices)) {
    x = indices[i,1]
    y = indices[i,2]
    values = c(values, z[(x-1):(x+1),(y-1):(y+1)])
  } ## i
  if (any(values < z[x0,y0])) {
    zz = 0
  } else {
    zz = 1
  }
  list(minima = zz, indices = indices)
}

imregionalmin = function (x) {
  
  is.minima = function(x, i, j) {
    mask = c(1,2,3,4,6,7,8,9)
    xx = x[(i-1):(i+1),(j-1):(j+1)]
    if (any(xx[mask] < xx[5])) {
      z = 0
    } else if (all(xx[mask] > xx[5])) {
      z = 1
    } else {
      zs = is.connected.minima(x, i, j)
      z = zs$minima
      zz[zs$indices] <<- zs$minima
    }
    return(z)
  }
  
  ## Expand grid to take care of edge cases
  xx = matrix(max(x) + 1, nrow(x) + 2, ncol(x) + 2)
  xx[2:(nrow(x)+1),2:(ncol(x)+1)] = x

  ## Find minima
  zz = matrix(-1, nrow(x) + 2 , ncol(x) + 2)
  zz[,1] = zz[1,] = zz[nrow(zz),] = zz[,ncol(zz)] = 0
  for (i in 1:nrow(x)) {
    ii = i + 1
    for (j in 1:ncol(x)) {
      #print(c(i,j))
      jj = j + 1
      if (zz[ii,jj] > -1)
        next
      zz[ii,jj] = is.minima(xx, ii, jj)
    } ## j
  } ## i
  
  ## Return result
  return(zz[2:(nrow(x)+1),2:(ncol(x)+1)])
  
}

## Find regional maxima
imregionalmax = function (x) {
  
  imregionalmin(-x)
  
}

## Bounding box check
point.in.bbox = function(x, y, bbox) {
  
  (bbox[1] <= x & x <= bbox[3]) & (bbox[2] <= y & y <= bbox[4])
  
}
