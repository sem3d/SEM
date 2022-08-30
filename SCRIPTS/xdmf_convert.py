import sys
import os
import os.path as osp
import numpy as np
from xml.etree import ElementTree as ET
from xml.etree import ElementInclude as EI
import h5py

srcres=sys.argv[1]
dstdir=sys.argv[2]
start=int(sys.argv[3])
stop=int(sys.argv[4])
step=int(sys.argv[5])

STEPS=range(start,stop,step)

FIELDS=["veloc"]
XFIELDS=["Veloc"]


def read_xdmf(srcres, xfields, steps):
    f=file(osp.join(srcres,"results_clean.xmf"))
    root = ET.fromstring(f.read())
    EI.include(root)
    
    domain = root.getchildren()[0]
    
    # The top spatial collection
    top_grid = domain.getchildren()[0]
    spatialgrids = top_grid.getchildren()
    
    newtop=ET.Element("Grid", GridType= 'Collection',
                      CollectionType='Temporal', Name='root')
    
    temporal_grids = [ sp.getchildren() for sp in spatialgrids ]
    
    
    common = {}
    
    ts = 1
    for sg in zip(*temporal_grids):
        if sg[0].tag == "DataItem":
            name = sg[0].attrib["Name"]
            common[name] = sg
            continue
        time = sg[0].find("Time")
        spgrid = newtop.makeelement("Grid", { "GridType" : 'Collection', "CollectionType" : 'Spatial', "Name" : 'timestep.%d' % ts })
        spgrid.append(time)
        for i, subgrid in enumerate(sg):
            ng = spgrid.makeelement("Grid", subgrid.attrib)
            for item in subgrid.getchildren():
                if item.tag == "Time":
                    continue
                if item.tag == "Attribute":
                    name = item.attrib["Name"]
                    if name in common:
                        item.getchildren()[0] = common[name][i]
                    if name not in xfields:
                        continue
                ng.append(item)
            spgrid.append(ng)
        if ts in steps:
            newtop.append(spgrid)
        ts = ts+1
    
    domain.getchildren()[0] = newtop
    return root

def get_results(srcres):
    res = []
    geom = []
    for fname in os.listdir(srcres):
        if fname.startswith("Rsem"):
            res.append(fname)
        if fname.startswith("geometry"):
            geom.append(fname)
    nprocs = len(geom)
    return res, nprocs

def get_steps(srcres, steps):
    res, nprocs = get_results(srcres)
    r2 = []
    for k in steps:
        rdir = "Rsem%04d" %k
        if rdir in res:
            r2.append(rdir)
    return r2, nprocs


def extract_fields(srcres, dstdir, step, nprocs):
    rdir = "Rsem%04d" % step
    fullrdir = osp.join(dstdir,rdir)
    if not os.access(fullrdir, os.F_OK):
        os.makedirs(fullrdir)
    print "Extracting:", fullrdir
    for k in range(nprocs):
        srcfield = osp.join(srcres,rdir,"sem_field.%04d.h5" % k)
        dstfield = osp.join(dstdir,rdir,"sem_field.%04d.h5" % k)
        print srcfield, "->", dstfield
        f = h5py.File(srcfield,"r")
        g = h5py.File(dstfield,"w")
        for nm in FIELDS:
            g.create_dataset(nm, shape=f[nm].shape, dtype=np.float32, compression=6)
            if f[nm].shape[0]!=0:
                data = f[nm][...]
                g[nm][...] = data
        f.close()
        g.close()
    print "done"



if __name__ == "__main__":
    res, nprocs = get_steps(srcres, STEPS)
    root = read_xdmf(srcres, XFIELDS, STEPS)
    output = osp.join(dstdir, "res2.xmf")
    file(output,"w").write(ET.tostring(root))
    for step in STEPS:
        extract_fields(srcres, dstdir, step, nprocs)



