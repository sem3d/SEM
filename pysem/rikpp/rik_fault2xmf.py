# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Write xmf file for fault points
"""
u'''Required modules'''
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import h5py
from rik_pp_lib import *
import xdmf_help as xdmf

u'''General informations'''
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2018, CentraleSupélec (MSSMat UMR CNRS 8579)"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"


def write_fault_xdmf(filename,Name_of_h5,c,tEnd,dtm,nel,nnd):
    fid = open(filename, 'w')
    
    # Header for xml file
    xdmf.write_header(fid)
    t=0
    if c=='moment':
        pfx = 'mom'
    elif c=='slip':
        pfx = 'slp'
    elif c=='srate':
        pfx = 'sra'
    while t <= tEnd:
        attr = [{'{}_{:>d}'.format(pfx,t):(c,'Node','Scalar',(nnd,))}]
        xdmf.write_geometry(fid,Name_of_h5,t*dtm,nel=nel,nnd=nnd,tag='fault',gtp='XYZ')
        xdmf.write_attributes(fid,Name_of_h5,attr)
        fid.write('''</Grid>\n''')
        t = t + 1
    xdmf.write_close(fid)
    fid.close()

class cfault2xmf(object):
    def __init__(self,opt):
       self.opt=opt 
    def convert(self):
        mf = opj(self.opt['wkd'],self.opt['mf'])
        kf = opj(self.opt['wkd'],self.opt['kf'])

        dic = {'moment':mf,'slip':kf,'srate':kf}
        Name_of_h5 = opj(self.opt['wkd'],mf.split('/')[-1].split('.')[0].split('_')[0])
        
        with h5py.File(mf,"r+") as mf, h5py.File(kf,"r+") as kf:
            # Fault grid
            nodes = mf['Nodes'][...]
            elems = mf['Elements'][...]
            
            nnd = nodes.shape[0]
            nel = elems.shape[0]
            
            dt = kf.attrs["dt"] 
            Nt = kf.attrs["Nt"]
            # Total time steps
            tEnd = Nt*dt
    
            for c in self.opt['ct']:
                filename = '{0}_{1}.xmf'.format(Name_of_h5,c)
                write_fault_xdmf(filename,dic[c],c,Nt,dt,nel,nnd)
    
if __name__=='__main__':
    opt = start_rik()
    fx = cfault2xmf(opt)
    fx.convert()
