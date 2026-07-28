# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Main file to read h5 mesh files
"""
#===============================================================================
# Required modules
#===============================================================================
# General modules
import h5py
import numpy as np
from . import mesh_shapes as mshp
#===============================================================================
# General informations
#===============================================================================
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2018, CentraleSupélec (MSSMat UMR CNRS 8579)"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

#===============================================================================
# SHAPE DEFINITION
#===============================================================================
#sph = mshp.sphere(1550.,0.,-1550.,200.)
sph = mshp.ellipsoid(1000,0.,-800.,400.,100.,200.)
#===============================================================================
# PROPERTIES
#===============================================================================
props = ['Nu','Vs','Rho']
names = {'Nu':'/workdir/gattif/Sl-Fl-3D-Het-Obj/mat-het-ell/h5/Mat_1_Nu.h5',
         'Vs':'/workdir/gattif/Sl-Fl-3D-Het-Obj/mat-het-ell/h5/Mat_1_Vs.h5',
        'Rho':'/workdir/gattif/Sl-Fl-3D-Het-Obj/mat-het-ell/h5/Mat_1_Density.h5'}
#values = {'Kappa':100.*10**9,'Mu':65.*10**9,'Rho':7850.}
values = {'Nu':0.48,'Vs':8000.,'Rho':7850.}
#===============================================================================
# OVERWRITE HDF5 FILE
#===============================================================================
for prop in props:
    name = names[prop]
    value= values[prop]
    fmesh = h5py.File(name,"r+")
    xMinGlob = fmesh.attrs["xMinGlob"]
    xMaxGlob = fmesh.attrs["xMaxGlob"]
    xStep    = fmesh.attrs["xStep"]
    xn = np.arange(xMinGlob[0],xMaxGlob[0]+xStep[0],xStep[0],float)
    yn = np.arange(xMinGlob[1],xMaxGlob[1]+xStep[1],xStep[1],float)
    zn = np.arange(xMinGlob[2],xMaxGlob[2]+xStep[2],xStep[2],float)
    xv,yv,zv = np.meshgrid(xn,yn,zn,sparse=False,indexing='xy')
    mat = fmesh["samples"][...]
    mat = np.rollaxis(mat,1,0)
    mat = np.rollaxis(mat,2,1)
    zone = sph.iswithin(xv,yv,zv)
    mat[zone] = value
    mat = np.rollaxis(mat,2,1)
    mat = np.rollaxis(mat,1,0)
    del fmesh['samples']
    dset = fmesh.create_dataset('samples', data=mat)
    fmesh.close()
