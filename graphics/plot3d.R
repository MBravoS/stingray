# NB: first complie "routines.R" before running this file

# stingray filename
file = '/Users/do/Dropbox/Code/Fortran/stingray/stingray/example/mock.hdf5'

# load data
gal = h5read(file,'galaxies',bit64conversion='bit64')
group = h5read(file,'groups',bit64conversion='bit64')
tile = h5read(file,'tiles')
shell = h5read(file,'shells')
para = h5read(file,'parameters')
n = length(gal$mag)

# make cartesian coordinates
gal$x = sph2car(gal$ra,gal$dec,gal$dc) # [Mpc/h]

# 3D plot
rgl.closeall()
rgl.sphgrid(radius=as.numeric(max(shell$dc_max)*para$box_side),radaxis = FALSE, col.long = 'black', col.lat = 'black')
rgl.tiling(tile,shell,para$box_side)
#rgl.shells(shell,para$box_side)
color = rev(rainbow(1000,start=0,end=2/3))
f = as.integer(gal$dc)
col = color[(f-min(f))/(max(f)-min(f))*999+1]
points3d(gal$x, col = col)
