# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Write xmf file for fault points
"""
#===============================================================================
# Required modules
#===============================================================================
# General modules
from __future__ import division
import warnings
warnings.filterwarnings("ignore")
import os.path as osp
import argparse
import numpy as np
import h5py
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


def write_xmf(filename,Name_of_h5,contour,tEnd,waittime,nOutput,Nx):
    f = open(filename, 'w')
    
    # Header for xml file
    f.write('''<?xml version="1.0" ?>
    <!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
    <Xdmf Version="2.0">
    <Domain>
    <Grid Name="Box" GridType="Collection" CollectionType="Temporal">
    ''')
    
    # loop over the attributes name written using time
    t = 0
    #dx_plane = xgrid[0,1]-xgrid[0,0]
    #dy_plane = ygrid[1,0]-ygrid[0,0]
    #dz_plane = depth[0,1]-depth[0,0]
    #print('dx-dy-dz plane')
    #print(dx_plane,dy_plane,dz_plane)
    #dx_rot = np.dot(MatMesh,np.array([dx_plane,dy_plane,dz_plane]))
    #print('dx-dy-dz rot')
    #print(dx_rot[0],dx_rot[1],dx_rot[2])

    dx_rot = np.array([25.791207949560572, 106.75152144602042, 89.1005859375])
    dts = ['x','y','z']
    frameN = 0 # For time sequence 
    while t <= tEnd :
        t = t + 1; 
        if( np.mod(t, nOutput) == 0 and t > waittime):
    
            # Naming datasets 
            dataSetName1 = '%s.h5:/%s'%(Name_of_h5,contour)
            #dataSetName2 = '%s.h5:/V_%.8d'%(Name_of_h5,t)
    
            # at individual time write the time independent Box grid. is it overdoing?
            f.write('''
            <!-- time step -->
            <Grid Name="Box %d" GridType="Uniform"> # 
            <Topology TopologyType="2DCoRectMesh" Dimensions="%d %d"/>
            <Geometry GeometryType="ORIGIN_DXDYDZ">
               <DataItem DataType="Float" Dimensions="3" Format="XML">
               %15.5f
               %15.5f
               %15.5f
               </DataItem>
               <DataItem DataType="Float" Dimensions="3" Format="XML">
               %15.5f
               %15.5f
               %15.5f
               </DataItem>
            </Geometry>
            <Time Value="%d" />
            '''%(frameN, Nx[0], Nx[1],dx_rot[0],dx_rot[1],dx_rot[2],
                dx_rot[0],dx_rot[1],dx_rot[2],frameN))
    
            # First Attribute
            f.write('''\n
            <Attribute Name="%s" AttributeType="Scalar" Center="Node">
            <DataItem Dimensions="%d %d" NumberType="Float" Precision="4"
            Format="HDF">%s
            </DataItem>
            </Attribute>
            '''%(contour,Nx[0],Nx[1],dataSetName1))
    
            ## Second Attribute
            #f.write('''\n
            #<Attribute Name="N" AttributeType="Vector" Center="Node">
            #<DataItem Dimensions="%d %d %d 3" NumberType="Float" Precision="4"
            #Format="HDF"> %s
            #</DataItem>
            #</Attribute>
            #</Grid>\n'''%(Nx, Ny, Nz, dataSetName2))
            frameN +=1
    
    # End the xmf file
    f.write('''
       </Grid>
    </Domain>
    </Xdmf>
    ''')


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # Parse database folder
    parser.add_argument('-m','--moment',
                        type  = str, 
                        default =  'napa2014_moment.h5',  
                        help = 'Moment rate tensor file')
    parser.add_argument('-k','--kine', 
                        type  = str, 
                        default =  'napa2014_kine.h5', 
                        help = 'Fault parameters file')
    parser.add_argument('-c','--contours',
                        nargs='*',
                        default = ['moment'], 
                        help = 'Fault parameters file')

    parser.add_argument('-s','--snapshots',
                        type = float,
                        default = 10., 
                        help = 'Time interval between snapshots')
    args = parser.parse_args()
    #
    mfn = args.moment
    kfn = args.kine
    Name_of_h5 = mfn.split('/')[-1].split('.')[0]
    # Contours
    contours = args.contours
    print('Contours')
    print(contours)
    # Data dumping step
    nOutput = args.snapshots
    # Not taking initial points
    waittime = 0
    
    with h5py.File(mfn,"r+") as mf, h5py.File(kfn,"r+") as kf:
        # Fault grid
        xg = kf["x"][...] 
        yg = kf["y"][...]
        zg = kf["z"][...]
        Nx = xg.shape 
        Ny = yg.shape
        Nz = zg.shape
        print("Fult Grid:")
        print("X = {}; Y = {}; Z = {}".format(Nx,Ny,Nz))
#[TODO] ADD STRIKE/DIP/SLIP TO COMPUTE ROT MATRIX
        ## Les angles 
        #strike =
        #dip    =
        #rake   =
        ## Les angles en radian
        #aS    = np.pi/180.0* strike
        #aR    = np.pi/180.0* rake
        #aD    = np.pi/180.0* dip
        ## Rotation matrix
        #MatMesh = rp.get_rotation_tensor(aS, aD)
        # Time steps
        dt = kf.attrs["dt"] 
        Nt = kf.attrs["Nt"]
        print("Time specifics:")
        print("dt = {}; Nt = {}".format(dt,Nt))
        # Total time steps
        tEnd = Nt*dt

        for c in contours:
            filename = '{0}_{1}.xmf'.format(Name_of_h5,c)
            write_xmf(filename,Name_of_h5,c,tEnd,waittime,nOutput,Nx)
