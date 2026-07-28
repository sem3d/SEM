# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Script read SEM3D snapshots with MPI-based algorithm

    Ex.1 : Parse Velocity snapshot files with n MPI cores
        
        mpirun --np n python3 parse_h5_snapshots.py @@wkd /path/to/sem3d/res @@var veloc

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
import debugpy
import pyvista as pv

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
       'Jac':'static','Kappa':'static','Lamb':'static','Material':'static',
       'Mu':'static','Nodes':'static','Proc':'static','displ':'dynamic','veloc':'dynamic','Accel':'dynamic',
       'eps_dev_xx':'dynamic','eps_dev_xy':'dynamic','eps_dev_xz':'dynamic',
       'eps_dev_yy':'dynamic','eps_dev_yz':'dynamic','eps_dev_zz':'dynamic',
       'sig_dev_xx':'dynamic','sig_dev_xy':'dynamic','sig_dev_xz':'dynamic',
       'sig_dev_yy':'dynamic','sig_dev_yz':'dynamic','sig_dev_zz':'dynamic',
       'eps_vol':'dynamic','press_elem':'dynamic','press_gll':'dynamic'}

sup = {'Mass':'node','Dens':'node','Dom':'node','Elements':'element',
       'Jac':'node','Kappa':'node','Lamb':'node','Material':'element',
       'Mu':'node','Nodes':'node','Proc':'element',
       'press_gll':'node','displ':'node','veloc':'node','Accel':'node',
       'eps_dev_xx':'element','eps_dev_xy':'element','eps_dev_xz':'element',
       'eps_dev_yy':'element','eps_dev_yz':'element','eps_dev_zz':'element',
       'sig_dev_xx':'element','sig_dev_xy':'element','sig_dev_xz':'element',
       'sig_dev_yy':'element','sig_dev_yz':'element','sig_dev_zz':'element',
       'eps_vol':'element','press_elem':'element'}

# def intersect_indices(x, y):
#     u_x, u_idx_x = np.unique(x, return_index=True)
#     u_y, u_idx_y = np.unique(y, return_index=True)
#     i_xy = np.intersect1d(u_x, u_y, assume_unique=True)
#     i_idx_x = u_idx_x[np.in1d(u_x, i_xy, assume_unique=True)]
#     i_idx_y = u_idx_y[np.in1d(u_y, i_xy, assume_unique=True)]
#     return i_idx_x, i_idx_y
def GetNonUniqueIndexes(UniqueArray,Array):
    sorted_keys = np.argsort(UniqueArray)
    indexes = sorted_keys[np.searchsorted(UniqueArray, Array, sorter=sorted_keys)]
    return indexes

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
        self.snapfile['geo'] = glob.glob(osj(self.wkd,'geometry*.h5'))
        self.snapfile['res'] = glob.glob(osj(self.wkd,'Rsem*'))
        self.snapfile['nc']  = len(self.snapfile['geo'])
        if self.end_time == -1: 
            self.end_time = len(self.snapfile['res'])
        self.nt = self.end_time - self.begin_time + 1
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

    def NodeCoordinates2Hash(self,NodeCoords):
        return np.array([hashlib.md5(i.tobytes()).digest() for i in NodeCoords], dtype="S16")
    
    def ParseSEM3DSnapshots(self):
        
        self.GlobalReNumbering()
        
        # if self.flag['static']:
        for g in self.snapfile['np']:
            with hf.File(osj(self.wkd,'geometry{:>04d}.h5'.format(g)),'r') as h5f:
                for v in self.flag['static']:
                    if self.dset[v].size == 0:
                        self.dset[v]=h5f[v][...]
                    else:
                        self.dset[v]=np.append(self.dset[v],h5f[v][...],axis=0)
                    
        for v in self.flag['static']:
            if sup[v]=="node":
                self.dset[v] = self.dset[v][self.LocalElementConnectivityOriginal]
                self.dset[v] = self.dset[v][self.LocalNodeUniqueHashIndex]
                self.dset[v] = self.dset[v][self.Local2UniqueLocalIndexOnRank]
        
        
        for v in self.flag['dynamic']:
            for g in self.snapfile['np']:
                with hf.File(osj(self.wkd,'Rsem{:>04d}/sem_field.{:>04d}.h5'.format(self.begin_time,g)),'r') as h5f:
                    if self.dset[v].size == 0:
                        self.dset[v]=h5f[v][...]
                    else:
                        self.dset[v]=np.append(self.dset[v],h5f[v][...],axis=0)
            
            if len(self.dset[v].shape)==2:
                self.dset[v]=self.dset[v].reshape((*self.dset[v].shape,1))
            elif len(self.dset[v].shape)==1:
                self.dset[v]=self.dset[v].reshape((*self.dset[v].shape,1,1))
            if sup[v]=="node":
                self.dset[v] = self.dset[v][self.LocalElementConnectivityOriginal]
                self.dset[v] = self.dset[v][self.LocalNodeUniqueHashIndex]
                self.dset[v] = self.dset[v][self.Local2UniqueLocalIndexOnRank]
            shp = self.dset[v].shape
            
            for t in range(self.begin_time+1, self.end_time+1):
                dtmp = np.array([])
                for g in self.snapfile['np']:
                    with hf.File(osj(self.wkd,'Rsem{:>04d}/sem_field.{:>04d}.h5'.format(t,g)),'r') as h5f:
                        if dtmp.size == 0:
                            dtmp=h5f[v][...]
                        else:
                            dtmp=np.append(dtmp,h5f[v][...],axis=0)
                if sup[v]=="node":
                    dtmp = dtmp[self.LocalElementConnectivityOriginal]
                    dtmp = dtmp[self.LocalNodeUniqueHashIndex]
                    dtmp = dtmp[self.Local2UniqueLocalIndexOnRank]
                self.dset[v]=np.append(self.dset[v], dtmp.reshape(*shp),axis=2)    
    
    def GlobalReNumbering(self):
        
        LocalNodeCoordinates = None
        LocalElementConnectivity = None
        NodeCount = 0
        self.ElementCount = 0
        # Loop over geometry.XXXX.h5 files (static data) on process
        # Might be more than one file per process!
        for g in self.snapfile['np']:
            with hf.File(osj(self.wkd,'geometry{:>04d}.h5'.format(g)),'r') as h5f:
                # Local node coordinates truncated at 6 decimals
                OnFileNodeCoordinates = h5f["Nodes"][...].astype(np.float64).round(decimals=6)
                if LocalNodeCoordinates is None:
                    LocalNodeCoordinates = OnFileNodeCoordinates
                else:
                    LocalNodeCoordinates = np.append(LocalNodeCoordinates,
                                                     OnFileNodeCoordinates,
                                                     axis=0)
                
                # Element connectivity on geometry file
                OnFileElementConnectivity = h5f["Elements"][...].astype(np.int64)+NodeCount
                
                if LocalElementConnectivity is None:
                    LocalElementConnectivity = OnFileElementConnectivity
                else:
                    LocalElementConnectivity = np.append(LocalElementConnectivity,
                                                         OnFileElementConnectivity,
                                                         axis=0)
                NodeCount += OnFileNodeCoordinates.shape[0]
                self.ElementCount += OnFileElementConnectivity.shape[0]

        num_elements = LocalElementConnectivity.shape[0]
        cells = []
        for conn in LocalElementConnectivity:
            cells.append(8)
            cells.extend(conn)
        cell_types = np.full(num_elements, 12, dtype=np.uint8)
        cells = np.array(cells,dtype=np.int64)
        hex_mesh = pv.UnstructuredGrid(cells,cell_types, LocalNodeCoordinates)
        hex_mesh.save(f"multi_hexahedral_mesh_{self.rank}.vtk")
        
        LocalNodeHash = self.NodeCoordinates2Hash(LocalNodeCoordinates)

        LocalNodeUniqueHash, self.LocalNodeUniqueHashIndex, LocalNodeUniqueHashInverse = np.unique(LocalNodeHash[LocalElementConnectivity.flatten()],
                                                                                                   return_index=True,
                                                                                                   return_inverse=True,
                                                                                                   axis=0)
        self.LocalElementConnectivityOriginal = LocalElementConnectivity.flatten().copy()
        LocalNodeCoordinates = LocalNodeCoordinates[LocalElementConnectivity.flatten(),:][self.LocalNodeUniqueHashIndex,:]
        self.LocalNodeCount = LocalNodeUniqueHash.size
        LocalElementConnectivity = LocalNodeUniqueHashInverse.reshape(-1,8)
        
        
        # num_elements = LocalElementConnectivity.shape[0]
        # cells = []
        # for conn in LocalElementConnectivity:
        #     cells.append(8)
        #     cells.extend(conn)
        # cell_types = np.full(num_elements, 12, dtype=np.uint8)
        # cells = np.array(cells,dtype=np.int64)
        # hex_mesh = pv.UnstructuredGrid(cells,cell_types, LocalNodeCoordinates)
        # hex_mesh.save(f"multi_hexahedral_mesh_unique_{self.rank}.vtk")
        
        GlobalNodeHash = self.comm.allgather(LocalNodeUniqueHash)
        RankPerGlobalNode = np.concatenate([r*np.ones((GlobalNodeHash[r].size,)) for r in range(self.size)],
                                           axis=0).astype(np.int64)

        GlobalNodeHash, GlobalUniqueIndex = np.unique(np.concatenate(GlobalNodeHash),
                                                      axis=0,
                                                      return_index=True)
        self.GlobalNumberofNodes = GlobalNodeHash.size
        RankPerGlobalNode = RankPerGlobalNode[GlobalUniqueIndex]
        OnRankIndexPerGlobalNode = np.argwhere(RankPerGlobalNode==self.rank).flatten()
        
        _, Local2UniqueLocalIndexOnRank, Global2UniqueLocalIndexOnRank = np.intersect1d(LocalNodeUniqueHash,
                                                                                        GlobalNodeHash,
                                                                                        assume_unique=True,
                                                                                        return_indices=True)
        _, index, _  = np.intersect1d(Global2UniqueLocalIndexOnRank,
                                     OnRankIndexPerGlobalNode,
                                     assume_unique=True,
                                     return_indices=True)
  
        self.dset["ElementsGlobal"] = Global2UniqueLocalIndexOnRank[LocalElementConnectivity]
        self.dset["NodesGlobal"] = LocalNodeCoordinates[Local2UniqueLocalIndexOnRank[index]]
        self.Local2UniqueLocalIndexOnRank = Local2UniqueLocalIndexOnRank[index]
        self.Global2UniqueLocalIndexOnRank = Global2UniqueLocalIndexOnRank[index]
        self.GlobalNodeCount = self.dset["NodesGlobal"].shape[0]
        
        ### CHECKING
        # CheckElements = self.comm.allgather(self.dset["ElementsGlobal"].flatten())
        # CheckNodes = self.comm.allgather(self.dset["NodesGlobal"].flatten())
        # CheckIndices = self.comm.allgather(self.Global2UniqueLocalIndexOnRank)
        # CheckElements = np.concatenate(CheckElements,axis=0).reshape(-1,8)
        # CheckNodes1 = np.concatenate(CheckNodes,axis=0).reshape(-1,3)
        # CheckNodes = np.empty_like(CheckNodes1)
        # CheckNodes[np.concatenate(CheckIndices)]=CheckNodes1
        
        
        # num_elements = CheckElements.shape[0]
        # cells = []
        # for conn in CheckElements:
        #     cells.append(8)
        #     cells.extend(conn)
        # cell_types = np.full(num_elements, 12, dtype=np.uint8)
        # cells = np.array(cells,dtype=np.int64)
        # hex_mesh = pv.UnstructuredGrid(cells,cell_types, CheckNodes)
        # hex_mesh.save(f"multi_hexahedral_mesh_unique_check_{self.rank}.vtk")
        #######
    

def ParseCL():
    """
        Parse command line flags
    """
    parser = argparse.ArgumentParser(prefix_chars='@')
    parser.add_argument('@@wkd',type=str,default='./res',help="Path to res directory")
    parser.add_argument('@@var',type=str,nargs='+',default=['Mass','Jac','Mu','Lamb',
                                                            'Elements',
                                                            'Dom','displ',
                                                            'eps_vol','eps_dev_xx',
                                                            'eps_dev_yy','eps_dev_zz',
                                                            'eps_dev_xy','eps_dev_yz',
                                                            'eps_dev_xz'],
                        help="Select snapshot")
    parser.add_argument('@@begin_time','@b', type=int, default=1, help="Initial time step to parse")
    parser.add_argument('@@end_time','@e', type=int, default=-1, help="Final time step to parse")
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
    snp.ParseSEM3DSnapshots()

    return snp

def main():
    """
    Start parallel process
    """
    MPI.Init()
    comm = MPI.COMM_WORLD               # Get communicator
    size = MPI.COMM_WORLD.Get_size()    # Get size of communicator
    rank = MPI.COMM_WORLD.Get_rank()    # Get the current rank
    
    # comm.Barrier()
    # port = 5678 + rank
    # debugpy.log_to("./")
    # debugpy.listen(("localhost", port))
    # debugpy.wait_for_client()
    # print(f"Rank {rank} waiting for debugger attach at port {port}")
    
    snp = GetSnapshots(comm,size,rank)

    # Gather global elements
    NumberofElements = np.array(comm.allgather(snp.ElementCount*8))
    
    if rank == 0:
        recvbuf = np.empty(sum(NumberofElements),dtype=np.int64)
    else:
        recvbuf = None
    
    comm.Gatherv(sendbuf=snp.dset["ElementsGlobal"].flatten(), 
                    recvbuf=[recvbuf, NumberofElements,
                            [sum(NumberofElements[:i]) for i in range(len(NumberofElements))],
                            MPI.LONG_LONG],
                    root=0)
    
    if rank == 0:
        ElementsGlob = recvbuf.astype(np.int64).reshape((-1,8))
    else:
        ElementsGlob = None
        
    
    
    # Gather global nodes
    NumberofNodeCoordinates = np.array(comm.allgather(snp.GlobalNodeCount*3))

    if rank == 0:
        recvbuf = np.empty(sum(NumberofNodeCoordinates),dtype=np.float64)
    else:
        recvbuf = None
    comm.Gatherv(sendbuf=snp.dset["NodesGlobal"].astype(np.float64).flatten(), 
                 recvbuf=[recvbuf, NumberofNodeCoordinates,
                          [sum(NumberofNodeCoordinates[:i]) for i in range(len(NumberofNodeCoordinates))],
                          MPI.DOUBLE],
                 root=0)

    if rank == 0:
        NodesTemp = recvbuf.astype(np.float64).reshape((-1,3))
    else:
        NodesTemp = None
        
    # Gather node global inindexes
    NumberofNodes = np.array(comm.allgather(snp.GlobalNodeCount))
    if rank == 0:
        recvbuf = np.empty(sum(NumberofNodes),dtype=np.float64)
    else:
        recvbuf = None
    comm.Gatherv(sendbuf=snp.Global2UniqueLocalIndexOnRank.astype(np.float64), 
                 recvbuf=[recvbuf, NumberofNodes,
                          [sum(NumberofNodes[:i]) for i in range(len(NumberofNodes))],
                          MPI.DOUBLE],
                 root=0)


    if rank == 0:
        NodesIndex = recvbuf.astype(np.int64)
        NodesGlobal = np.empty_like(NodesTemp)
        NodesGlobal[NodesIndex] = NodesTemp
    else:
        NodesIndex = None
    
    NumberofVectorComponents = NumberofNodeCoordinates*snp.nt
    if rank == 0:
        recvbuf = np.empty(sum(NumberofVectorComponents),dtype=np.float64)
    else:
        recvbuf = None
    comm.Gatherv(sendbuf=snp.dset["displ"].astype(np.float64).flatten(), 
                 recvbuf=[recvbuf, NumberofVectorComponents,
                          [sum(NumberofVectorComponents[:i]) for i in range(len(NumberofVectorComponents))],
                          MPI.DOUBLE],
                 root=0)

    if rank == 0:
        DisplTemp = recvbuf.astype(np.float64).reshape((-1,3,snp.nt))
        DisplGlobal = np.empty_like(DisplTemp)
        DisplGlobal[NodesIndex] = DisplTemp
    else:
        DisplGlobal = None
        
    
    if rank == 0:
        import pyvista as pv
        # Convert to VTK-compatible format
        num_elements = ElementsGlob.shape[0]
        cells = []
        for conn in ElementsGlob:
            cells.append(8)  # First value: number of nodes per hexahedron
            cells.extend(conn)  # Add node indices
        cells = np.array(cells, dtype=np.int64)  # Convert to NumPy array

        # Cell types (all VTK_HEXAHEDRON)
        cell_types = np.full(num_elements, 12, dtype=np.uint8)
        # Create the Unstructured Grid
        hex_mesh = pv.UnstructuredGrid(cells,cell_types, NodesGlobal)
        # hex_mesh.save("multi_hexahedral_mesh.vtk")
        # hex_mesh["Time"] = np.arange(snp.begin_time, snp.end_time+1)
        for j in range(DisplGlobal.shape[-1]):
            t = snp.begin_time + j
            hex_mesh["displ"] = DisplGlobal[:,:,j].reshape((-1,3))
            hex_mesh.save(f"multi_hexahedral_mesh_{t}.vtu")
        
    MPI.Finalize()

if __name__=="__main__":
    main()