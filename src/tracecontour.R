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
  
  ## Expand domain to allow for objects on top and bottom rows
  xx = matrix(0, nx, ny + 2)
  xx[,2:(ny + 1)] = x
  ny = ny + 2
  
  ## Find first meridian with no points in it
  ## Need to start from leftmost point in object, this avoids finding the
  ## right side of a wrapped object which is stored on the left.
  i0 = 0
  total = 1
  while(total > 0 & i0 <= nx) {
    i0 = i0 + 1
    total = sum(x[i0,])
  }
  i0 = (i0 - 1) %% nx + 1
  
  ## Find first filled point, scanning from meridians from bottom to top and
  ## left to right
  c = c(1,1)
  skip = FALSE
  for (i in i0:nx) {
    for (j in 1:ny) {
      if (x[i,j] == 1) {
        s = c(i,j)
        skip = TRUE
        break
      } else {
        c = c(i,j)
      }
    } ## j
    if (skip)
      break
  } ## i
  
  ## Trace contour using Moore-Neighbor tracing
  p = s
  B = p
  c0 = c
  cm1 = c
  c = move.clockwise(c, p)
  cm1x = (cm1[1] - 1) %% nx + 1
  cm1y = cm1[2]
  cx = (c[1] - 1) %% nx + 1
  cy = c[2]
  while(cx != s[1] | cy != s[2] | cm1x != c0[1] | cm1y != c0[2]) {
    ## If
    if (x[cx,cy] == 1) {
      B = cbind(B,c)
      p = c
      c = cm1
      c = move.clockwise(c, p)
    } else {
      cm1 = c
      cm1x = (cm1[1] - 1) %% nx + 1
      cm1y = cm1[2]
      c = move.clockwise(c, p)
    }
    cx = (c[1] - 1) %% nx + 1
    cy = c[2]
  } ## c != s
  
  ## Remove added latitude and return
  B[2,] = B[2,] - 1
  t(B)
  
}
