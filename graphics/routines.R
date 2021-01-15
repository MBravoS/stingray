# libraries
library(rhdf5)
library(rgl)
library(sphereplot)
library(bit64)
library(magicaxis)

# routines
rgl.new = function(windowSize = 700, bg = 'white') {
  open3d(windowRect=c(0,0,windowSize,windowSize))
  bg3d(col=bg)
}

rgl.closeall = function() {
  while(rgl.cur()[[1]]>0) {rgl.close()}
}

rgl.hold = function() par3d(skipRedraw=TRUE)

rgl.done = function() par3d(skipRedraw=FALSE)

rgl.tiling = function(tile,shell,zoom=1,col='grey') {
  par3d(skipRedraw=TRUE) # stop drawing
  for (itile in seq(length(tile$tile_id))) {
    ishell = tile$shell[itile]
    sxx = shell$transformation$rotation_xx[ishell]
    sxy = shell$transformation$rotation_xy[ishell]
    sxz = shell$transformation$rotation_xz[ishell]
    syx = shell$transformation$rotation_yx[ishell]
    syy = shell$transformation$rotation_yy[ishell]
    syz = shell$transformation$rotation_yz[ishell]
    szx = shell$transformation$rotation_zx[ishell]
    szy = shell$transformation$rotation_zy[ishell]
    szz = shell$transformation$rotation_zz[ishell]
    rotationmatrix = rbind(c(sxx,sxy,sxz),c(syx,syy,syz),c(szx,szy,szz))
    rgl.cube(tile$center_x[itile],tile$center_y[itile],tile$center_z[itile],rotationmatrix,as.numeric(zoom),col=col)
  }
  par3d(skipRedraw=FALSE) # draw now
}

rgl.shells = function(shell,zoom=1) {
  par3d(skipRedraw=TRUE) # stop drawing
  if (shell$dc_min[1]>0) rgl.sphere(0,0,0, radius = shell$dc_min[1]*as.numeric(zoom))
  for (i in seq_along(shell$shell_id)) {
    rgl.sphere(0,0,0, radius = as.numeric(shell$dc_max[i])*as.numeric(zoom),col='#000000',alpha=0.05)
  }
  par3d(skipRedraw=FALSE) # draw now
}

rgl.cube = function(x,y,z,rotationmatrix=diag(c(1,1,1)),zoom=1,col='grey') {
  rot.lines3d(x+c(0.5,0.5,-0.5,-0.5,0.5),y+c(0.5,-0.5,-0.5,0.5,0.5),z+c(0.5,0.5,0.5,0.5,0.5),rotationmatrix,zoom,col=col)
  rot.lines3d(x+c(0.5,0.5,-0.5,-0.5,0.5),y+c(0.5,-0.5,-0.5,0.5,0.5),z-c(0.5,0.5,0.5,0.5,0.5),rotationmatrix,zoom,col=col)
  rot.lines3d(x+c(0.5,0.5),y+c(0.5,0.5),z+c(-0.5,0.5),rotationmatrix,zoom,col=col)
  rot.lines3d(x+c(0.5,0.5),y-c(0.5,0.5),z+c(-0.5,0.5),rotationmatrix,zoom,col=col)
  rot.lines3d(x-c(0.5,0.5),y+c(0.5,0.5),z+c(-0.5,0.5),rotationmatrix,zoom,col=col)
  rot.lines3d(x-c(0.5,0.5),y-c(0.5,0.5),z+c(-0.5,0.5),rotationmatrix,zoom,col=col)
}

rot.lines3d = function(x,y,z,rot,zoom,...) {
  x = rot%*%rbind(x,y,z)*zoom
  lines3d(t(x),...)
}

rgl.cone = function(para, nalpha = 20, nr = 20) {
  
  par3d(skipRedraw=TRUE) # stop drawing
  
  # initialize variables
  angle = para$angle
  axis = c(para$axis.x,para$axis.y,para$axis.z)
  dcmin = para$dc_min/para$L
  dcmax = para$dc_max/para$L
  L = para$L
  dmin = para$dc_min/para$L
  dmax = para$dc_max/para$L
  axis = axis/sqrt(sum(axis^2))
  
  # central line
  lines3d(axis[1]*c(dmin,dmax),axis[2]*c(dmin,dmax),axis[3]*c(dmin,dmax),col = 'grey')
  
  # cone
  ex = c(1.23,9.843,-2.34827)
  e1 = vectProd(axis,ex)
  e1 = e1/vectNorm(e1)
  e2 = vectProd(axis,e1)
  e2 = e2/vectNorm(e2)
  dalpha = 2*pi/nalpha
  p = array(NA,c(4,3))
  v = array(NA,c(nr+1,3))
  for (alpha in seq(0,2*pi,length=nalpha+1)) {
    for (ir in seq(1,nr+1)) {
      v[ir,] = axis*cos(angle*(ir-1)/nr)+(cos(alpha)*e1+sin(alpha)*e2)*sin(angle*(ir-1)/nr)
    }
    if (alpha>0) {
      
      # cone surface
      p[1,] = dmin*v[nr+1,]
      p[2,] = dmax*v[nr+1,]
      p[3,] = dmax*v.old[nr+1,]
      p[4,] = dmin*v.old[nr+1,]
      polygon3d(p,col='grey',alpha=0.2)
      
      # cone end
      for (ir in seq(1,nr)) {
        p[1,] = dmax*v[ir,]
        p[2,] = dmax*v[ir+1,]
        p[3,] = dmax*v.old[ir+1,]
        p[4,] = dmax*v.old[ir,]
        polygon3d(p,col='grey',alpha=0.2)
      }
      
      # cone start
      if (dmin>0) {
        for (ir in seq(1,nr)) {
          p[1,] = dmin*v[ir,]
          p[2,] = dmin*v[ir+1,]
          p[3,] = dmin*v.old[ir+1,]
          p[4,] = dmin*v.old[ir,]
          polygon3d(p,col='grey',alpha=0.2)
        }
      }
      
    }
    v.old = v
  }
  
  par3d(skipRedraw=FALSE) # draw now
}

rgl.sphere <- function (x, y=NULL, z=NULL, ng=50, radius = 1, color="white", add=F, ...) {
  lat <- matrix(seq(90, -90, len = ng)*pi/180, ng, ng, byrow = TRUE)
  long <- matrix(seq(-180, 180, len = ng)*pi/180, ng, ng)
  
  vertex  <- rgl:::rgl.vertex(x, y, z)
  nvertex <- rgl:::rgl.nvertex(vertex)
  radius  <- rbind(vertex, rgl:::rgl.attr(radius, nvertex))[4,]
  color  <- rbind(vertex, rgl:::rgl.attr(color, nvertex))[4,]
  
  for(i in 1:nvertex) {
    x <- vertex[1,i] + radius[i]*cos(lat)*cos(long)
    y <- vertex[2,i] + radius[i]*cos(lat)*sin(long)
    z <- vertex[3,i] + radius[i]*sin(lat)
    persp3d(x, y, z, specular="white", add=TRUE, color=color[i], ...)
  }
}

vectProd <- function(x,y){
  # cross-product of two 3-vectors
  return(c(x[2]*y[3]-x[3]*y[2],x[3]*y[1]-x[1]*y[3],x[1]*y[2]-x[2]*y[1]))
}

vectNorm <- function(x){
  # norm of an n-vector
  return(sqrt(sum(x*x)))
}