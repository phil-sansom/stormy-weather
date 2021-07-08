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
  
  nx = nrow(x)
  ny = ncol(x)
  c = c(1,1)
  skip = FALSE
  for (i in 1:nx) {
    for (j in 1:ny) {
      if (buffer[i,j] == 1) {
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
    if (buffer[cx,cy] == 1) {
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
  
}
