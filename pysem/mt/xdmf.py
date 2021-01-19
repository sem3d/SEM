"""
Utility function to create
Xdmf file from a set of meshes descriptions
"""
import xml.etree.ElementTree as ET
import h5py
from mt.mesh_files import Points, Mesh, ELEMENT_TYPES
import os.path as osp



TYPE_DICT = {
    "float32" : ("Float", 4),
    "float64" : ("Float", 8),
    "int8"    : ("Int"  , 1),
    "int16"   : ("Int"  , 2),
    "int32"   : ("Int"  , 4),
    "int64"   : ("Int"  , 8),
    "uint8"   : ("UInt" , 1),
    "uint16"  : ("UInt" , 2),
    "uint32"  : ("UInt" , 4),
    "uint64"  : ("UInt" , 8),
}

ATTRIBUTE_TYPE = {
    1 : "Scalar",
    3 : "Vector",
    6 : "Tensor6",
    9 : "Tensor"
}

XDMF_HEADER ="""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">
"""

def get_type(dtype):
    return TYPE_DICT[dtype.name]

class XdmfWriter(object):
    def __init__(self):
        self.xdmf = ET.Element("Xdmf", attrib = {"Version":"2.0",
                                                 "xmlns:xi":"http://www.w3.org/2001/XInclude"})
        self.root = ET.ElementTree(self.xdmf)
        self.dom = ET.SubElement(self.xdmf, "Domain")
        self.root.text = "\n"
        self.xdmf.text = "\n"
        self.dom.text = "\n"
        self.xdmf.tail = "\n"

    def write(self, f):
        with open(f,"w") as fid:
            fid.write(XDMF_HEADER)
        with open(f,"wb") as fid:
            self.root.write(fid, xml_declaration=False)

    def uniform_grid(self, parent, name):
        return ET.SubElement(parent, "Grid", attrib={"Name":name, "GridType":"Uniform" })

    def spatial_grid(self, parent, name, temporal=False):
        if temporal:
            ctype = "Temporal"
        else:
            ctype = "Spatial"
        return ET.SubElement(parent, "Grid",
                             attrib = { "Name" : name,
                                        "GridType" : "Collection",
                                        "CollectionType" : ctype })

    def data_item(self, parent, fname, path, dset):
        dims = " ".join(["%s" % d for d in dset.shape])
        dtype, prec = get_type(dset.dtype)
        item = ET.SubElement(parent, "DataItem",
                          attrib={ "Format" : "HDF",
                                   "Dimensions" : dims,
                                   "NumberType" : dtype,
                                   "Precision" : str(prec) })
        item.text = "%s:%s" % (fname, path)
        item.tail = "\n"
        return item

    def attribute(self, parent, name, dset, attr_type="Cell"):
        dims = " ".join(["%s" % d for d in dset.shape])
        nc = 1
        if len(dset.shape)==2:
            nc = dset.shape[1]
        item = ET.SubElement(parent, "Attribute",
                             attrib={ "Name" : name,
                                      "Center" : attr_type,
                                      "AttributeType" : ATTRIBUTE_TYPE.get(nc, "Matrix"),
                                      "Dimensions" : dims
                                  }
                         )
        item.text="\n"
        item.tail="\n"
        return item

    def data_item_ref(self, parent, ref):
        pass

    def geometry(self, parent, typ="XYZ"):
        item = ET.SubElement(parent, "Geometry",
                             attrib = { "Type" : typ } )
        item.text="\n"
        item.tail="\n"
        return item

    def topology(self, parent, cells):
        topo_type = cells.topology_type()
        item = ET.SubElement(parent, "Topology",
                             attrib = { "Type" : topo_type,
                                        "NumberOfElements" : str(cells.nelems()),
                                    })
        item.text = "\n"
        item.tail = "\n"
        return item

    def structured_topology(self, parent, topo_type, dims):
        item = ET.SubElement(parent, "Topology",
                             attrib = { "Type" : topo_type,
                                        "Dimensions" : " ".join("%s" %d for d in dims),
                             })
        item.text = "\n"
        item.tail = "\n"
        return item

    def time(self, grid, val):
        item = ET.SubElement(grid, "Time",
                             attrib = {"Value" : str(val)})
        item.tail="\n"
        return item
        

def write_xdmf(fname, glob_nodes, groups, node_fields, temporal=False):
    w = XdmfWriter()
    parent = w.dom
    nodes_name = "/Nodes"
    if len(groups)>1:
        # Create a spatial collection
        parent = w.spatial_grid(parent, "main", temporal)
    for (elems, grp) in groups:
        fname = grp.file.filename
        path = grp.name
        grid = w.uniform_grid(parent, "main")
        if "Time" in grp.attrs:
            w.time(grid, grp.attrs["Time"])
        if "Nodes" in grp:
            nodes = grp["Nodes"]
            nodes_name = osp.join(path,"Nodes")
        else:
            nodes = glob_nodes
            nodes_name = "/Nodes"
        elems.write_geometry(w, fname, grid, grp, nodes_name, nodes)
        elems.write_topology(w, fname, grid, grp)
        # Cell/Node data
        for k in grp.keys():
            attr_type = None
            if k == "Nodes":
                continue
            if k=="Mat":
                attr_type = "Cell"
                typ = "Cell"
                name = "Mat"
            if k.startswith("Cell"):
                attr_type = "Cell"
                typ, name = k.split("_",1)
            if k.startswith("Node"):
                attr_type = "Node"
                typ, name = k.split("_",1)
            if attr_type is None:
                continue
            data = grp[k]
            attr = w.attribute(grid, name, data, attr_type=attr_type)
            itm = w.data_item(attr, fname, osp.join(grp.name, k), data)
        # Global nodes
        for dname, dset in node_fields:
            typ, name = dname.split("_",1)
            attr = w.attribute(grid, name, dset, attr_type="Node")
            itm = w.data_item(attr, fname, dname, dset)
    
    w.write(fname+".xmf")


def create_xdmf_structure(fname, temporal=False):
    f = h5py.File(fname,"r")
    groups = []
    node_fields = []
    nodes = None
    if "Nodes" in f:
        nodes = f["Nodes"]
    for grpname in f.keys():
        grp = f[grpname]
        if isinstance(grp, h5py.Dataset):
            if grpname!="Nodes" and grpname.startswith("Node"):
                node_fields.append((grpname, grp))
                continue
        if not isinstance(grp, h5py.Group):
            continue
        for typ in ELEMENT_TYPES:
            if typ.match_group(grp):
                break
        else:
            print("Group {} doesn't contain cell information".format(grpname))
            # Does not contain any cell description leave it...
            continue
        elems = typ.from_group(grp)
        groups.append((elems, grp))
    write_xdmf(fname, nodes, groups, node_fields, temporal)
