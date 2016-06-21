import numpy as np
import os.path as osp

class KlassName(object):
    def __get__(self, obj, cls=None):
        return cls.__name__

class BaseElement(object):
    """Base class for the description of a mesh element, contains
    two class members :
    name : the name of the element type, appears as a dataset name of elements of that type
    nnodes : the number of nodes for element of that types, -1 if variable number
    in which case, a second field "Indices" should be present within the HDF5 file"""
    def get_name(self):
        return self.__class__.__name__

    name = KlassName()
    nnodes = None
    xdmf_topology_type = "unknown"

    def __init__(self, elements, indices=None):
        self.elements = elements
        self.indices = indices
        if self.nnodes>0:
            assert self.elements.shape[1] == self.nnodes
        elif self.nnodes==-1:
            assert indices is not None

    def nelems(self):
        return self.elements.shape[0]

    @classmethod
    def match_group(cls, grp):
        """Tells whether an HDF5 group contains data to match these element types"""
        return cls.name in grp.keys()
    
    @classmethod
    def from_group(cls, grp):
        elems = grp[cls.name]
        if cls.nnodes == -1:
            indices = grp["Indices"]
        else:
            indices = None
        return cls(elems, indices)

    def topology_type(self):
        """Return a type string for xml format"""
        return self.xdmf_topology_type

    def write_topology(self, w, fname, grid, grp):
        topo = w.topology(grid, self)
        data = w.data_item(topo, fname, osp.join(grp.name, self.name), self.elements)

    def write_geometry(self, w, fname, grid, grp, nodes_name, nodes):
        geom = w.geometry(grid, typ="XYZ")
        data = w.data_item(geom, fname, nodes_name, nodes)


class Seg2(BaseElement):
    nnodes = 2
    xdmf_topology_type = "Polyline"

class Tri3(BaseElement):
    nnodes = 3
    xdmf_topology_type = "Triangle"

class Quad4(BaseElement):
    nnodes = 4
    xdmf_topology_type = "Quadrilateral"

class Quad8(BaseElement):
    nnodes = 8
    xdmf_topology_type = "Quadrilateral_8"

class Poly2d(BaseElement):
    nnodes = -1

class Tetra4(BaseElement):
    nnodes = 4
    xdmf_topology_type = "Tetrahedron"

class Hexa8(BaseElement):
    nnodes = 8
    xdmf_topology_type = "Hexahedron"

class Hexa27(BaseElement):
    nnodes = 27
    xdmf_topology_type = "Hexahedron_27"

class Poly3d(BaseElement):
    nnodes = -1

class StructGrid2D(BaseElement):
    xdmf_topology_type = "2DSMesh"

    @classmethod
    def match_group(cls, grp):
        """Tells whether an HDF5 group contains data to match these element types"""
        if "X" not in grp or "Y" not in grp:
            return False
        if len(grp["X"].shape)!=2 or len(grp["Y"].shape)!=2:
            return False
        return grp["X"].shape == grp["Y"].shape

    @classmethod
    def from_group(cls, grp):
        vx = grp["X"]
        vy = grp["Y"]
        vz = grp["Z"]
        elems = (vx,vy,vz)
        indices = None
        return cls(elems, indices)

    def nelems(self):
        gx = self.elements[0]
        return (gx.shape[0]-1)*(gx.shape[1]-1)

    def write_topology(self, w, fname, grid, grp):
        elshape = [ d for d in self.elements[0].shape ]
        topo = w.structured_topology(grid, self.xdmf_topology_type, elshape)
        #data = w.data_item(topo, fname, osp.join(grp.name, self.name), self.elements)

    def write_geometry(self, w, fname, grid, grp, nodes_name, nodes):
        geom = w.geometry(grid, typ="X_Y_Z")
        namex = osp.join(grp.name, "X")
        data = w.data_item(geom, fname, namex, self.elements[0])
        namex = osp.join(grp.name, "Y")
        data = w.data_item(geom, fname, namex, self.elements[1])
        namex = osp.join(grp.name, "Z")
        data = w.data_item(geom, fname, namex, self.elements[2])

class StructGrid3D(BaseElement):
    @classmethod
    def match_group(cls, grp):
        """Tells whether an HDF5 group contains data to match these element types"""
        if "X" not in grp or "Y" not in grp or "Z" not in grp:
            return False
        if len(grp["X"].shape)!=3 or len(grp["Y"].shape)!=3 or len(grp["Z"].shape)!=3:
            return False
        return (grp["X"].shape == grp["Y"].shape) and (grp["X"].shape == grp["Z"].shape)

    @classmethod
    def from_group(cls, grp):
        vx = grp["X"]
        vy = grp["Y"]
        vz = grp["Z"]
        elems = (vx,vy,vz)
        indices = None
        return cls(elems, indices)


ELEMENT_TYPES = [ Seg2, Tri3, Quad4, Quad8, Poly2d, Tetra4, Hexa8, Hexa27, Poly3d, StructGrid2D, StructGrid3D ]

class Points(object):
    """Manages allocation of unique 3d points"""
    def __init__(self, PREC=10000.):
        self.prec = PREC
        self.points = {}

    def count(self):
        return len(self.points)

    def add(self, x, y, z=0.):
        """This functions maintains a dictionary of points coordinate
        in order to associate a unique number to each point"""
        prec = self.prec
        pt = (int(x*prec), int(y*prec), int(z*prec))
        num = len(self.points)
        return self.points.setdefault(pt, num)

    def add_array(self, x, y, z=None):
        if z is None:
            z = np.zeros_like(x)
        prec = self.prec
        aix = np.array(x*prec, int)
        aiy = np.array(y*prec, int)
        aiz = np.array(z*prec, int)
        pnum = np.zeros_like(x, int)
        for n,(i,j,k) in enumerate(zip(aix.flat, aiy.flat, aiz.flat)):
            pt = (int(i),int(j),int(k))
            num = len(self.points)
            pnum.flat[n] = self.points.setdefault(pt, num)
        return pnum

    def get_nodes(self):
        NP = len(self.points)
        prec = float(self.prec)
        P = np.zeros((NP,3), float)
        for (ix, iy, iz), n in self.points.items():
            P[n,0] = ix
            P[n,1] = iy
            P[n,2] = iz
        P /= prec
        return P


class Mesh(object):
    def __init__(self, PREC=10000.):
        self.points = Points(PREC)
        self.quad4 = []
        self.quad8 = []
        self.tri3 = []
        self.hexa8 = []
        self.hexa27 = []
        self.tetra4 = []

    def add_pt(self, x, y, z=0.):
        return self.points.add(x,y,z)

    # TODO accessors ? add_quad4, add_quad9, ...
    def add_quad4_grid(self, X, Y, Z=None):
        """Add a structured grid as unstructured quad4"""
        ni, nj = X.shape
        if Z is None:
            Z = np.zeros_like(X)
        for i in range(ni-1):
            for j in range(nj-1):
                n0 = self.points.add(X[i  ,j  ], Y[i  ,j  ], Z[i  ,j  ])
                n1 = self.points.add(X[i  ,j+1], Y[i  ,j+1], Z[i  ,j+1])
                n2 = self.points.add(X[i+1,j+1], Y[i+1,j+1], Z[i+1,j+1])
                n3 = self.points.add(X[i+1,j  ], Y[i+1,j  ], Z[i+1,j  ])
                self.quad4.append( (n0,n1,n2,n3) )

    def add_quad8_grid(self, X, Y, Z=None):
        """Add a structured grid as unstructured quad4"""
        ni, nj = X.shape
        if Z is None:
            Z = np.zeros_like(X)
        for i in range(0, ni-1, 2):
            for j in range(0, nj-1, 2):
                n0 = self.points.add(X[i  ,j  ], Y[i  ,j  ], Z[i  ,j  ])
                n1 = self.points.add(X[i  ,j+2], Y[i  ,j+2], Z[i  ,j+2])
                n2 = self.points.add(X[i+2,j+2], Y[i+2,j+2], Z[i+2,j+2])
                n3 = self.points.add(X[i+2,j  ], Y[i+2,j  ], Z[i+2,j  ])
                n4 = self.points.add(X[i  ,j+1], Y[i  ,j+1], Z[i  ,j+1])
                n5 = self.points.add(X[i+1,j+2], Y[i+1,j+2], Z[i+1,j+2])
                n6 = self.points.add(X[i+2,j+1], Y[i+2,j+1], Z[i+2,j+1])
                n7 = self.points.add(X[i+1,j  ], Y[i+1,j  ], Z[i+1,j  ])
                self.quad8.append( (n0,n1,n2,n3,n4,n5,n6,n7) )

    def add_hexa8_grid(self, X, Y, Z):
        """Add a structured grid as unstructured quad4"""
        ni, nj, nk = X.shape
        for i in range(ni-1):
            for j in range(nj-1):
                for k in range(nk-1):
                    n0 = self.points.add(X[i  ,j  ,k  ], Y[i  ,j  ,k  ], Z[i  ,j  ,k  ])
                    n1 = self.points.add(X[i  ,j+1,k  ], Y[i  ,j+1,k  ], Z[i  ,j+1,k  ])
                    n2 = self.points.add(X[i+1,j+1,k  ], Y[i+1,j+1,k  ], Z[i+1,j+1,k  ])
                    n3 = self.points.add(X[i+1,j  ,k  ], Y[i+1,j  ,k  ], Z[i+1,j  ,k  ])
                    n4 = self.points.add(X[i  ,j  ,k+1], Y[i  ,j  ,k+1], Z[i  ,j  ,k+1])
                    n5 = self.points.add(X[i  ,j+1,k+1], Y[i  ,j+1,k+1], Z[i  ,j+1,k+1])
                    n6 = self.points.add(X[i+1,j+1,k+1], Y[i+1,j+1,k+1], Z[i+1,j+1,k+1])
                    n7 = self.points.add(X[i+1,j  ,k+1], Y[i+1,j  ,k+1], Z[i+1,j  ,k+1])
                    self.hexa8.append( (n0,n1,n2,n3,n4,n5,n6,n7) )

    def add_hexa8_gridn(self, N):
        """Add a structured grid as unstructured quad4"""
        ni, nj, nk = N.shape
        for i in range(ni-1):
            for j in range(nj-1):
                for k in range(nk-1):
                    n0 = N[i  ,j  ,k  ]
                    n1 = N[i  ,j+1,k  ]
                    n2 = N[i+1,j+1,k  ]
                    n3 = N[i+1,j  ,k  ]
                    n4 = N[i  ,j  ,k+1]
                    n5 = N[i  ,j+1,k+1]
                    n6 = N[i+1,j+1,k+1]
                    n7 = N[i+1,j  ,k+1]
                    self.hexa8.append( (n0,n1,n2,n3,n4,n5,n6,n7) )

    
    def write_nodes(self, f, nodes=None):
        if "Nodes" in f:
            del f["Nodes"]
        if nodes is None:
            nodes = self.points.get_nodes()
        f["Nodes"] = nodes

    def write_elements(self, f, name, data, dsetname):
        grp = f.get(name, None)
        if grp is None:
            grp = f.create_group(name)
        data = np.array(data, np.uint64)
        grp.create_dataset(dsetname, dtype=np.uint64, data=data,compression=6,shuffle=True)
        
    def write_quad4(self, f, name, **kwargs):
        self.write_elements(f, name, self.quad4, "Quad4")

    def write_quad8(self, f, name, **kwargs):
        self.write_elements(f, name, self.quad8, "Quad8")

    def write_tri3(self, f, name, **kwargs):
        self.write_elements(f, name, self.tri3, "Tri3")

    def write_hexa8(self, f, name, **kwargs):
        self.write_elements(f, name, self.hexa8, "Hexa8")

    def write_hexa27(self, f, name, **kwargs):
        self.write_elements(f, name, self.hexa27, "Hexa27")

    def write_tetra4(self, f, name, **kwargs):
        self.write_elements(f, name, self.tetra4, "Tetra4")

    def get_cell(self, f, name):
        grp = f[name]
        for typ in ELEMENT_TYPES:
            if typ.name in grp:
                return typ.from_group(grp)
        return None

    def write_cell_data(self, f, name, dsetname, data):
        grp = f[name]
        cells = self.get_cell(f, name)
        assert cells is not None
        assert data.shape[0] == cells.nelems()
        if len(data.shape)==1:
            nc = 1
        else:
            nc = data.shape[1]
        # XXX ? compression ?
        grp.create_dataset(dsetname, data=data)

    def write_node_data(self, f, dsetname, data):
        if len(data.shape)==1:
            nc = 1
        else:
            nc = data.shape[1]
        # XXX ? compression ?
        f.create_dataset(dsetname, data=data)

