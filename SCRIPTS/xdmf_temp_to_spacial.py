import sys

from xml.etree import ElementTree as ET
from xml.etree import ElementInclude as EI
f=file(sys.argv[1])
root = ET.fromstring(f.read())
EI.include(root)

domain = root.getchildren()[0]

# The top spatial collection
top_grid = domain.getchildren()[0]
spatialgrids = top_grid.getchildren()

newtop=ET.Element("Grid", GridType= 'Collection', CollectionType='Temporal', Name='root')

temporal_grids = [ sp.getchildren() for sp in spatialgrids ]


common = {}

ts = 0
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
            ng.append(item)
        spgrid.append(ng)
    newtop.append(spgrid)
    ts = ts+1

domain.getchildren()[0] = newtop

print ET.tostring(root)

