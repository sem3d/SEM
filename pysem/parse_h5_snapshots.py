# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Script to compute PML properties for SEM3D

    Ex.1 : Compute amplitude Ax knowning PML length (300. m):
        
        python3 compute_pml_length.py @@PML_length 300.
    
    Ex.2 : Compute PML_length knowing the amplitude Ax (10.)
        
        python3 compute_pml_length.py @@Ax 10.
"""
# Required modules
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI
import glob
import argparse
from os.path import join as osj
import numpy as np
import h5py as hf


# General informations
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2020, CentraleSupélec (MSSMat UMR CNRS 8579)"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

typ = {'Mass':'static','Dens':'static','Dom':'static','Elements':'static',
       'Jac':'static','Kappa':'static','Lamb':'static','Mass':'static','Material':'static',
       'Mu':'static','Nodes':'static','Proc':'static','displ':'dynamic','veloc':'dynamic','accel':'dynamic',
       'eps_dev_xx':'dynamic','eps_dev_xy':'dynamic','eps_dev_xz':'dynamic',
       'eps_dev_yy':'dynamic','eps_dev_yz':'dynamic','eps_dev_zz':'dynamic',
       'sig_dev_xx':'dynamic','sig_dev_xy':'dynamic','sig_dev_xz':'dynamic',
       'sig_dev_yy':'dynamic','sig_dev_yz':'dynamic','sig_dev_zz':'dynamic',
       'eps_vol':'dynamic','press_elem':'dynamic','press_gll':'dynamic'}

class snapshots(object):
    def __init__(self,**kwargs):
        self.__call__(**kwargs)
        
    def __call__(self,**kwargs):
        self.__dict__.update(**kwargs)
        self.setup()
        
    def setup(self):
        self.snapfile={}
        self.flag={}
        self.dset = {}
        self.dtmp = {}
        # self.comm = opt["comm"]
        # self.size = opt["size"]
        # self.rank = opt["rank"]
        self.snapfile['geo'] = glob.glob(osj(self.wkd,'geometry*.h5'))
        self.snapfile['res'] = glob.glob(osj(self.wkd,'Rsem*'))
        self.snapfile['nc']  = len(self.snapfile['geo'])
        self.snapfile['nt']  = len(self.snapfile['res'])
        (qc,rc) = divmod(self.snapfile['nc'],self.size)

        if self.rank<rc:
            self.snapfile['np'] = [self.rank*(qc+1)+q for q in range(qc+1)]
        else:
            self.snapfile['np'] = [rc*(qc+1)+(self.rank-rc)*qc+q for q in range(qc)]

        print("Rank {:d} - file range {}".format(self.rank,self.snapfile['np']))

        self.flag['static'] = []
        self.flag['dynamic']= []
        
        for v in self.var:
            if 'static' in typ[v]:
                if self.flag['static']:
                    self.flag['static'].append(v)
                else:
                    self.flag['static'] = [v]
            if 'dynamic' in typ[v]:
                if self.flag['dynamic']:
                    self.flag['dynamic'].append(v)
                else:
                    self.flag['dynamic'] = [v]
            self.dset[v] = np.array([])
    
    def parse(self):
        if self.flag['static']:
            for g in self.snapfile['np']:
                with hf.File(osj(self.wkd,'geometry{:>04d}.h5'.format(g)),'r') as fid:
                    for v in self.flag['static']:
                        if self.dset[v].size == 0:
                            self.dset[v]=fid[v][...]
                        else:
                            self.dset[v]=np.append(self.dset[v],fid[v][...],axis=0)
        if self.flag['dynamic']:
            for v in self.flag['dynamic']:
                for g in self.snapfile['np']:
                    with hf.File(osj(self.wkd,'Rsem0001/sem_field.{:>04d}.h5'.format(g)),'r') as fid:
                        if self.dset[v].size == 0:
                            self.dset[v]=fid[v][...]
                        else:
                            self.dset[v]=np.append(self.dset[v],fid[v][...],axis=0) 
                if len(self.dset[v].shape)==2:
                    self.dset[v]=self.dset[v].reshape((*self.dset[v].shape,1))
                elif len(self.dset[v].shape)==1:
                    self.dset[v]=self.dset[v].reshape((*self.dset[v].shape,1,1))
                shp = self.dset[v].shape
                for t in range(2,self.snapfile['nt']+1):
                    self.dtmp = np.array([])
                    for g in self.snapfile['np']:
                        with hf.File(osj(self.wkd,'Rsem{:>04d}/sem_field.{:>04d}.h5'.format(t,g)),'r') as fid:
                            if self.dtmp.size == 0:
                                self.dtmp=fid[v][...]
                            else:
                                self.dtmp=np.append(self.dtmp,fid[v][...],axis=0)
                    self.dset[v]=np.append(self.dset[v],self.dtmp.reshape(*shp),axis=2)

def ParseCL():
    """
        Parse command line flags
    """
    parser = argparse.ArgumentParser(prefix_chars='@')
    parser.add_argument('@@wkd',type=str,default='./res',help="Path to res directory")
    parser.add_argument('@@var',type=str,nargs='+',default=['Mass','Jac','Mu','Lamb',
                                                            'eps_dev_xz'],help="Select snapshot")
    opt = parser.parse_args().__dict__
    
    return opt

def GetSnapshots(comm,size,rank):

    # Parse Command Line
    opt = ParseCL()
    opt["comm"] = comm
    opt["size"] = size
    opt["rank"] = rank
    
    # Generate snapshot structure
    snp = snapshots(**opt)
    
    # Parse result snapshots 
    snp.parse()

    return snp

def main():
    """
    Start parallel process
    """
    MPI.Init()
    comm = MPI.COMM_WORLD               # Get communicator
    size = MPI.COMM_WORLD.Get_size()    # Get size of communicator
    rank = MPI.COMM_WORLD.Get_rank()    # Get the current rank
    hostname = MPI.Get_processor_name() # Get the hostname
    
    GetSnapshots(comm,size,rank)

    MPI.Finalize()

if __name__=="__main__":
    main()
    
#     eps_xx = snp.dset['eps_dev_xx']+snp.dset['eps_vol']/3.
#     eps_yy = snp.dset['eps_dev_yy']+snp.dset['eps_vol']/3.
#     eps_zz = snp.dset['eps_dev_zz']+snp.dset['eps_vol']/3.
#     eps_xy = snp.dset['eps_dev_xy']
#     eps_xz = snp.dset['eps_dev_xz']
#     eps_yz = snp.dset['eps_dev_yz']