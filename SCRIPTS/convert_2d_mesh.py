import sys
import h5py
import numpy as np

oldmesh = h5py.File(sys.argv[1])
newmesh = h5py.File(sys.argv[2],"w")

X = oldmesh["X"][...]
Y = oldmesh["Y"][...]
Z = oldmesh["Z"][...]
Q = oldmesh["Q"][...]
M1 = oldmesh["M1"][...]
M2 = oldmesh["M2"][...]

grp = newmesh.create_group("Sem2D")

nodes = np.zeros( (X.shape[0], 2), float )
nodes[:,0] = X
nodes[:,1] = Y

newmesh["Nodes"] = nodes

grp["Quad4"] = Q
grp["Mat"] = M1
grp["M2"] = M2


