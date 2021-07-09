move.clockwise = function(x, x0) {
  xmx0 = x - x0
  if (xmx0[1] == 0 & xmx0[2] == -1) {
    z = x0 + c(-1,-1)
  } else if (xmx0[1] == -1 & xmx0[2] == -1) {
    z = x0 + c(-1, 0)
  } else if (xmx0[1] == -1 & xmx0[2] ==  0) {
    z = x0 + c(-1, 1)
  } else if (xmx0[1] == -1 & xmx0[2] ==  1) {
    z = x0 + c( 0, 1)
  } else if (xmx0[1] ==  0 & xmx0[2] ==  1) {
    z = x0 + c( 1, 1)
  } else if (xmx0[1] ==  1 & xmx0[2] ==  1) {
    z = x0 + c( 1, 0)
  } else if (xmx0[1] ==  1 & xmx0[2] ==  0) {
    z = x0 + c( 1,-1)
  } else if (xmx0[1] ==  1 & xmx0[2] == -1) {
    z = x0 + c(0,-1)
  }
  # z[1] = (z[1] - 1) %% nx  + 1
  z
}

trace.contour = function(x) {
  
  ## Extract dimensions
  nx = nrow(x)
  ny = ncol(x)
  
  ## Expand domain to allow for points on top and bottom rows
  xx = matrix(0, nx, ny + 2)
  xx[,2:(ny + 1)] = x
  ny = ny + 2
  
  ## Find starting point, scanning from bottom to top and left to right
  i = 1
  search = TRUE
  while (i <= nx & search) {
    j = 2
    while (j <= ny - 1 & search) {
      if (xx[i,j] == 1) {
        start = c(i,j)
        search = FALSE
      }
      j = j + 1
    } ## j
    i = i + 1
  } ## i
  
  ## Check if there are more points to be found
  surround =   sum(xx[((start[1] - 1):(start[1] + 1) - 1) %% nx + 1,
                      (start[2] - 1):(start[2] + 1)])
  if (surround > 1) {
    
    ## Initialize storage
    boundary = matrix(NA, 0.5*nx*ny, 2)
    
    ## Trace contour using Moore-Neighbor tracing
    point = start                              ## Set pivot point to start point
    current = start                            ## Set search point to entry point
    current[2] = current[2] - 1
    previous = current                         ## Set fake previous search point
    current = move.clockwise(current, point)   ## Advance search point
    currentr = current
    currentr[1] = (currentr[1] - 1) %% nx + 1  ## Real coords of current point
    previousr = previous
    point1 = c(-1,-1)
    current1 = c(-1,-1)
    n = 0                                     ## Initialize counter
    while(any(currentr != point1) | any(previousr != current1)) {
      ## Check if point found
      if (xx[t(currentr)] == 1) {
        ## If point found store it, ...
        n = n + 1                      ## Increment counter
        point = current                ## Set pivot point to search point
        boundary[n,] = currentr        ## Record found point
        if (n == 1) {
          point1 = currentr
          current1 = previousr
        }
        current = previous             ## Set search point to entry point
        currentr = previousr
      } else {
        ## ... otherwise advance search
        previous  = current                       ## Set previous search point to current
        previousr = currentr
        current = move.clockwise(current, point)  ## Advance search point
        currentr = current
        currentr[1] = (currentr[1] - 1) %% nx + 1
      }
    } ## c != s
    
    ## Remove excess rows
    boundary = boundary[1:n,]
    
  } else {

    boundary = t(start)

  }
 
  ## Remove added latitude
  boundary[,2] = boundary[,2] - 1
     
  ## Return boundary
  return(boundary)
  
}
