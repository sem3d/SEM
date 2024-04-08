import h5py
import numpy as np

f=h5py.File("Mat_0.h5","w")
vp = f.create_group("Vp")
vs=f.create_group("Vs")
rho=f.create_group("Rho")

# Warning, dimensions are inversed due to fortran "ordering"
# and hdf5 library doing weird inversion stuff
nz,ny,nx=3,4,5
shp = nz,ny,nx

xvp = np.linspace(5000.,8000.,nx)
yvs = np.linspace(3000.,4000.,ny)
zrho = np.linspace(15000.,2200.,nz)
xmin = np.array( [-100.,-100.,-100] )
xmax = np.array( [500.,500.,500] )

# Vp
prop_vp=np.zeros( shp, float )
prop_vp[...] = xvp.reshape((1,1,-1))
vp["samples"] = prop_vp
vp.attrs["xMinGlob"] = xmin
vp.attrs["xMaxGlob"] = xmax

#Vs
prop_vs=np.zeros( shp, float )
prop_vs[...] = yvs.reshape((1,-1,1))
vs["samples"] = prop_vs
vs.attrs["xMinGlob"] = xmin
vs.attrs["xMaxGlob"] = xmax

#Rho
prop_rho=np.zeros( shp, float )
prop_rho[...] = zrho.reshape((-1,1,1))
rho["samples"] = prop_rho
rho.attrs["xMinGlob"] = xmin
rho.attrs["xMaxGlob"] = xmax

