import sys
import h5py
import numpy as np

oldmesh = h5py.File(sys.argv[1])
newmesh = h5py.File(sys.argv[2],"w")

nodes = oldmesh["nodes"][...]
Q = oldmesh["elements"][...]
M1 = oldmesh["material"][:,0]
M2 = oldmesh["material"][:,2]


grp = newmesh.create_group("Sem2D")

Nodes = np.zeros( (nodes.shape[0], 3), float)
print nodes.shape
Nodes[:,:2] = nodes
newmesh["Nodes"] = Nodes
grp["Quad4"] = Q
grp["Mat"] = M1
grp["M2"] = M2
grp["Cell_Mat"] = M1


