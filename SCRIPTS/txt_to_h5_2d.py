import sys
import h5py
import numpy as np

oldmesh = sys.argv[1]





newmesh = h5py.File(sys.argv[2],"w")


#X = oldmesh["X"][...]
#Y = oldmesh["Y"][...]
#Z = oldmesh["Z"][...]
#Q = oldmesh["Q"][...]
#M1 = oldmesh["M1"][...]
#M2 = oldmesh["M2"][...]

def read_intl(f):
    l = f.readline().split()
    return int(l[0])

def read_intls(f,n):
    l = f.readline().split()
    return [  int(v) for v in l[:n] ]

def read_floatl(f, n):
    l = f.readline().split()
    return [ float(v) for v in l[:n] ]

def read_mesh_txt(src):
    f = file(src)
    n = read_intl(f)
    nn = read_intl(f)
    nodes = np.zeros( (nn,3), float)
    for k in range(nn):
        nodes[k,:2] = read_floatl(f,2)

    nmat = read_intl(f)
    nel = read_intl(f)
    ncontrolnodes = read_intl(f)
    f.readline()
    Q = np.zeros( (nel,4), int)
    M1 = np.zeros( (nel,), int)
    M2 = np.zeros( (nel,), int)
    for k in range(nel):
        eldef = read_intls(f, 13)
        Q[k,:] = eldef[:4]
        M1[k] = eldef[12]
    return nodes, Q, M1, M2

nodes, Q, M1, M2 = read_mesh_txt(oldmesh)


grp = newmesh.create_group("Sem2D")


newmesh["Nodes"] = nodes
grp["Quad4"] = Q
grp["Mat"] = M1
grp["M2"] = M2
grp["Cell_Mat"] = M1


