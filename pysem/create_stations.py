# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Script to create stations' file for SEM3D

    Ex.1 : Create a [10 x 10 x 1] grid of receivers in a domain of 5m x 5m x 5m
        
        python3 create_stations.py -x -5. 5. -y -5. 5 -z 0. -s 10 10 1
"""
import argparse
import numpy as np

def grid(lims):
    xv = np.linspace(lims['xmin'],lims['xmax'],lims['nx'],dtype=np.float64)
    yv = np.linspace(lims['ymin'],lims['ymax'],lims['ny'],dtype=np.float64)
    zv = np.linspace(lims['zmin'],lims['zmax'],lims['nz'],dtype=np.float64)
    xg,yg,zg = np.meshgrid(xv,yv,zv,indexing='xy')
    return (xg,yg,zg)

def write_stations(xg,fn='stations.txt'):
    xg=tuple([x.squeeze().reshape(-1,1) for x in xg])
    np.savetxt(fn,np.concatenate(xg,axis=-1),fmt='%10.5f',delimiter=' ')

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlim",type=float,nargs="*",default=[-5.0,5.0],help="Limits of the box")
    parser.add_argument("--ylim",type=float,nargs="*",default=[-5.0,5.0],help="Limits of the box")
    parser.add_argument("--zlim",type=float,nargs="*",default=[0.0,0.0],help="Limits of the box")
    parser.add_argument("--step",type=float,nargs="*",default=[10,10,1],help="Numbers of points per direction")
    parser.add_argument("--indx",type=str,default='ij',help="Grid indexing xy|ij")
    parser.add_argument("--ofnm",type=str,default='stations.txt',help="Output filename")
    opt = parser.parse_args()

    assert len(opt.xlim)==2
    assert len(opt.ylim)==2
    assert len(opt.zlim)==2

    opt.xlim.sort() 
    opt.ylim.sort()
    opt.zlim.sort()
    
    lims = (('xmin',opt.xlim[0]),('xmax',opt.xlim[1]),
            ('ymin',opt.ylim[0]),('ymax',opt.ylim[1]),
            ('zmin',opt.zlim[0]),('zmax',opt.zlim[1]),
            ('nx',opt.step[0]),('ny',opt.step[1]),('nz',opt.step[2]))
    lims = dict(lims)
    print(f'Domain limits/discretization : \n X: {lims["xmin"]} m : {lims["xmax"]} m; nx={lims["nx"]} \n Y: {lims["ymin"]} m : {lims["ymax"]} m; ny={lims["ny"]} \n Z: {lims["zmin"]} m : {lims["zmax"]} m; nz={lims["nz"]}')
    
    # generate grid points
    grd = grid(lims)
    # write station file
    write_stations(grd,opt.ofnm)
    print("Station file {} has been created!".format(opt.ofnm))
