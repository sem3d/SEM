from __future__ import print_function
import sys


from xml.etree import ElementTree as ET
from xml.etree import ElementInclude as EI
f=file(sys.argv[1])

Nmax = -1
if len(sys.argv)>2:
    Nmax = int(sys.argv[2])

root = ET.fromstring(f.read())
EI.include(root)

domain = root.getchildren()[0]

# The top spatial collection
top_grid = domain.getchildren()[0]
spatialgrids = top_grid.getchildren()

newtop=ET.Element("Grid", GridType= 'Collection', CollectionType='Temporal', Name='root')

temporal_grids = [ sp.getchildren()[0].getchildren()[0] for sp in spatialgrids ]

common = {}

grids = {}

for temp_grid in temporal_grids:
    space_name = temp_grid.attrib["Name"]
    for mesh_grid in temp_grid.getchildren():
        mesh_name = mesh_grid.attrib["Name"]
        for i,el in enumerate(mesh_grid.getchildren()):
            if el.tag == "Time":
                break
        timeval = float(el.attrib["Value"].strip())
        del mesh_grid[i]
        _, gridlist = grids.setdefault(timeval,(el, []))
        gridlist.append(mesh_grid)

for ts, tm in enumerate(sorted(grids.keys())):
    timeel, gridlist = grids[tm]
    gridlist.sort(key=lambda g:g.attrib["Name"])
    spgrid = newtop.makeelement("Grid", { "GridType" : 'Collection', "CollectionType" : 'Spatial', "Name" : 'timestep.%d' % ts })
    spgrid.append(timeel)
    spgrid.tail="\n"
    for i,g in enumerate(gridlist):
        if i<3:
            continue
        spgrid.append(g)
        if Nmax>=0 and i>=Nmax:
            break
    newtop.append(spgrid)

## ts = 0
## for sg in temporal_grids:
##     if sg[0].tag == "DataItem":
##         name = sg[0].attrib["Name"]
##         common[name] = sg
##         continue
##     time = sg[0].find("Time")
##     if time is None:
##         toto
##     spgrid = newtop.makeelement("Grid", { "GridType" : 'Collection', "CollectionType" : 'Spatial', "Name" : 'timestep.%d' % ts })
##     spgrid.append(time)
##     for i, subgrid in enumerate(sg):
##         if Nmax>0 and i>=Nmax:
##             continue
##         ng = spgrid.makeelement("Grid", subgrid.attrib)
##         for item in subgrid.getchildren():
##             if item.tag == "Time":
##                 continue
##             if item.tag == "Attribute":
##                 name = item.attrib["Name"]
##                 if name in common:
##                     item.getchildren()[0] = common[name][i]
##             if item is None:
##                 continue
##             ng.append(item)
##         spgrid.append(ng)
##     newtop.append(spgrid)
##     ts = ts+1
## 
domain.getchildren()[0] = newtop

print(ET.tostring(root))

