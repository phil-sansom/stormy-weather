move.clockwise = function(x, x0, nlon) {
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
  # z[1] = (z[1] - 1) %% nlon  + 1
  z
}

trace.contour = function(x) {
  
  nx = nrow(x)
  ny = ncol(x)
  c = c(1,1)
  skip = FALSE
  for (ii in 1:nx) {
    for (jj in 1:ny) {
      if (buffer[ii,jj] == 1) {
        s = c(ii,jj)
        skip = TRUE
        break
      } else {
        c = c(ii,jj)
      }
    } ## jj
    if (skip)
      break
  } ## ii
  p = s
  B = p
  c0 = c
  cm1 = c
  c = move.clockwise(c, p, nx)
  cm1x = (cm1[1] - 1) %% nx + 1
  cm1y = cm1[2]
  cx = (c[1] - 1) %% nx + 1
  cy = c[2]
  while(cx != s[1] | cy != s[2] | cm1x != c0[1] | cm1y != c0[2]) {
    if (buffer[cx,cy] == 1) {
      B = cbind(B,c)
      p = c
      c = cm1
      c = move.clockwise(c, p, nx)
    } else {
      cm1 = c
      cm1x = (cm1[1] - 1) %% nx + 1
      cm1y = cm1[2]
      c = move.clockwise(c, p, nx)
    }
    cx = (c[1] - 1) %% nx + 1
    cy = c[2]
  } ## c != s
  B
  
}