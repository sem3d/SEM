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
import hashlib

# General informations
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2020, CentraleSupélec (MSSMat UMR CNRS 8579)"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

typ = {'Mass':'static','Dens':'static','Dom':'static','Elements':'static','ElementsGlob':'static',
       'Jac':'static','Kappa':'static','Lamb':'static','Mass':'static','Material':'static',
       'Mu':'static','Nodes':'static','Proc':'static','displ':'dynamic','veloc':'dynamic','accel':'dynamic',
       'eps_dev_xx':'dynamic','eps_dev_xy':'dynamic','eps_dev_xz':'dynamic',
       'eps_dev_yy':'dynamic','eps_dev_yz':'dynamic','eps_dev_zz':'dynamic',
       'sig_dev_xx':'dynamic','sig_dev_xy':'dynamic','sig_dev_xz':'dynamic',
       'sig_dev_yy':'dynamic','sig_dev_yz':'dynamic','sig_dev_zz':'dynamic',
       'eps_vol':'dynamic','press_elem':'dynamic','press_gll':'dynamic'}

# def intersect_indices(x, y):
#     u_x, u_idx_x = np.unique(x, return_index=True)
#     u_y, u_idx_y = np.unique(y, return_index=True)
#     i_xy = np.intersect1d(u_x, u_y, assume_unique=True)
#     i_idx_x = u_idx_x[np.in1d(u_x, i_xy, assume_unique=True)]
#     i_idx_y = u_idx_y[np.in1d(u_y, i_xy, assume_unique=True)]
#     return i_idx_x, i_idx_y

class SnapshotsSEM3D(object):
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

        self.flag['static']  = []
        self.flag['dynamic'] = []
        
        for v in self.var:
            if 'static' in typ[v] and v not in self.flag['static']:
                if self.flag['static']:
                    self.flag['static'].append(v)
                else:
                    self.flag['static'] = [v]
            if 'dynamic' in typ[v] and v not in self.flag['dynamic']:
                if self.flag['dynamic']:
                    self.flag['dynamic'].append(v)
                else:
                    self.flag['dynamic'] = [v]
            self.dset[v] = np.array([])

        self.global_renumbering()

    def global_renumbering(self):
        # count global node number
        lnnodes = []
        lnelems = 0
        lhashn = np.array([]).tobytes()
        for g in self.snapfile['np']:
            with hf.File(osj(self.wkd,'geometry{:>04d}.h5'.format(g)),'r') as h5f:
                # unique local node numbering
                lnodes = h5f["Nodes"][...].round(decimals=6)
                lnnodes.append(lnodes.shape[0])
                # node coordinate hash
                lhashn = np.append(lhashn,np.array([hashlib.md5(i).digest() for i in lnodes]))
                # Local Number of elements
                # lnelems += h5f["Elements"][...].shape[0]

        lhashn = np.unique(lhashn,axis=0)
        # global node hash
        self.ghashn = self.comm.allgather(lhashn)
        # global node hash
        self.ghashn = np.unique(np.concatenate(self.ghashn),axis=0)
        # globale number of nodes
        self.gnnodes = self.ghashn.size
        # global elements
        # self.lelems = np.empty((lnelems,8),dtype=np.int64)
        self.gnelems = self.comm.allreduce(lnelems,op=MPI.SUM)

        # Assemble global connectivity and unique nodal coordinates
        for g in self.snapfile['np']:
            with hf.File(osj(self.wkd,'geometry{:>04d}.h5'.format(g)),'r+') as h5f:
                # unique local node numbering
                lnodes = h5f["Nodes"][...].round(decimals=6)
                lconnect = h5f["Elements"][...]
                lnelems = lconnect.shape[0]
                hconnect = lnodes[lconnect].reshape((-1,3),order='C')
                hconnect = np.array([hashlib.md5(i).digest() for i in hconnect]).reshape(-1)
                
                for _ in range(8):
                    _,lnode, gnode = np.intersect1d(hconnect,self.ghashn,
                        assume_unique=False, return_indices=True)
                    lconnect[np.unravel_index(lnode,shape=(lnelems,8),order='C')] = gnode
                    hconnect[lnode] = None

                if "/ElementsGlob" not in h5f:
                    h5f.create_dataset("ElementsGlob", data=lconnect)
                else:
                    elemsglob = h5f["ElementsGlob"]
                    elemsglob[...] = lconnect
                h5f.close()

                    
    def parse(self):
        if self.flag['static']:
            for g in self.snapfile['np']:
                with hf.File(osj(self.wkd,'geometry{:>04d}.h5'.format(g)),'r') as h5f:
                    for v in self.flag['static']:
                        if self.dset[v].size == 0:
                            self.dset[v]=h5f[v][...]
                        else:
                            self.dset[v]=np.append(self.dset[v],h5f[v][...],axis=0)
        if self.flag['dynamic']:
            for v in self.flag['dynamic']:
                for g in self.snapfile['np']:
                    with hf.File(osj(self.wkd,'Rsem0001/sem_field.{:>04d}.h5'.format(g)),'r') as h5f:
                        if self.dset[v].size == 0:
                            self.dset[v]=h5f[v][...]
                        else:
                            self.dset[v]=np.append(self.dset[v],h5f[v][...],axis=0) 
                if len(self.dset[v].shape)==2:
                    self.dset[v]=self.dset[v].reshape((*self.dset[v].shape,1))
                elif len(self.dset[v].shape)==1:
                    self.dset[v]=self.dset[v].reshape((*self.dset[v].shape,1,1))
                shp = self.dset[v].shape
                for t in range(2,self.snapfile['nt']+1):
                    self.dtmp = np.array([])
                    for g in self.snapfile['np']:
                        with hf.File(osj(self.wkd,'Rsem{:>04d}/sem_field.{:>04d}.h5'.format(t,g)),'r') as h5f:
                            if self.dtmp.size == 0:
                                self.dtmp=h5f[v][...]
                            else:
                                self.dtmp=np.append(self.dtmp,h5f[v][...],axis=0)
                    self.dset[v]=np.append(self.dset[v],self.dtmp.reshape(*shp),axis=2)

def ParseCL():
    """
        Parse command line flags
    """
    parser = argparse.ArgumentParser(prefix_chars='@')
    parser.add_argument('@@wkd',type=str,default='./res',help="Path to res directory")
    parser.add_argument('@@var',type=str,nargs='+',default=['Mass','Jac','Mu','Lamb',
                                                            'Elements','ElementsGlob',
                                                            'Dom','displ',
                                                            'eps_vol','eps_dev_xx',
                                                            'eps_dev_yy','eps_dev_zz',
                                                            'eps_dev_xy','eps_dev_yz',
                                                            'eps_dev_xz'],
                        help="Select snapshot")
    opt = parser.parse_args().__dict__
    
    return opt

def GetSnapshots(comm,size,rank):

    # Parse Command Line
    opt = ParseCL()
    opt["comm"] = comm
    opt["size"] = size
    opt["rank"] = rank
    
    # Generate snapshot structure
    snp = SnapshotsSEM3D(**opt)
    
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
    
    snp = GetSnapshots(comm,size,rank)

    # unique_values,count = np.unique(snp.dset['ElementsGlob'],return_counts=True)
    # print(unique_values,count)
    MPI.Finalize()

if __name__=="__main__":
    main()
    
#     eps_xx = snp.dset['eps_dev_xx']+snp.dset['eps_vol']/3.
#     eps_yy = snp.dset['eps_dev_yy']+snp.dset['eps_vol']/3.
#     eps_zz = snp.dset['eps_dev_zz']+snp.dset['eps_vol']/3.
#     eps_xy = snp.dset['eps_dev_xy']
#     eps_xz = snp.dset['eps_dev_xz']
#     eps_yz = snp.dset['eps_dev_yz']