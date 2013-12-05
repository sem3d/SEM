# Ecriture d'un fichier mesh pour sem2d
import numpy as np
import h5py

class Mesh2D(object):
    def __init__(self, nodes):
        self.nodes = nodes
        self.elements = None
        self.materials = None
        self.el_edges = None
        self.edges = []
        self.rev_edges = {}

    def nmat(self):
        return np.unique(self.materials).size

    def nnodes(self):
        return self.nodes.shape[0]
    
    def nelements(self):
        return self.elements.shape[0]

    def nedges(self):
        return len(self.edges)

    def compute_connectivity(self, elements):
        self.elements = elements
        NE, _ = elements.shape
        self.el_edges = np.zeros( (NE,4), int)

        for i in range(NE):
            self.add_edge(i, 0, elements[i,0], elements[i,1])
            self.add_edge(i, 1, elements[i,1], elements[i,2])
            self.add_edge(i, 2, elements[i,3], elements[i,2])
            self.add_edge(i, 3, elements[i,0], elements[i,3])

        
    def add_edge(self, el, nf, n1, n2):
        if (n1,n2) in self.rev_edges:
            edge_num = self.rev_edges[ (n1,n2) ]
            edge_or = 0
            self.edges[edge_num][2] = (el,nf)
        elif (n2,n1) in self.rev_edges:
            edge_num = self.rev_edges[ (n2,n1) ]
            edge_or = 1
            self.edges[edge_num][2] = (el,nf)
        else:
            edge_num = len(self.rev_edges)
            edge_or = 0
            self.rev_edges[(n1,n2)] = edge_num
            self.edges.append( [(n1,n2), (el,nf), None] )

        self.el_edges[el,nf] = edge_num
        return edge_num, edge_or

def write_mesh(fname, mesh, proc):
    f = file(fname, "w")
    NN = mesh.nnodes()
    NE = mesh.nelements()
    elements = mesh.elements
    faces = mesh.el_edges
    materials = mesh.materials

    f.write("2\n") # Dimension
    f.write("%d\n" % NN) # Nb global nodes
    for i in range(NN):
        f.write("%f %f\n" % (mesh.nodes[i,0], mesh.nodes[i,1]) )
    
    f.write("%d\n" % mesh.nmat()) # Number of different materials
    f.write("%d\n" % NE)   # Number of local elements 
    f.write("4\n")   # number of nodes per element
    f.write("\n")
    for i in range(NE):
        f.write("%d %d %d %d %d %d %d %d %d %d %d %d %d\n" % (
            tuple(elements[i,:])+tuple(faces[i,:])+tuple(elements[i,:])+(materials[i],)))
                
    f.write("\n")
    f.write("%d\n" % mesh.nedges())
    f.write("\n")
    for i in range(mesh.nedges()):
        (n1, n2), (el1,nf1), xel2 = mesh.edges[i]
        if xel2 is not None:
            el2, nf2 = xel2
        else:
            el2, nf2 = -1, -1
        f.write("%d %d %d %d %d %d\n" % (el1, el2, nf1, nf2, n1, n2) )

    f.write("\n")
    f.write("%d\n" % NN)
    f.write("\n")
    for i in range(NN):
        f.write("%d\n" % i)
    f.write("\n\n1\n") # Number of partitions
    f.write("0\n") # number of communications


TMPL_H5="""<?xml version="1.0"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
<Domain>
<Grid CollectionType="Temporal" GridType="Collection">
  <DataItem Name="Mat" Format="HDF" Datatype="Int"  Dimensions="%(ne)s">%(fgeom)s:/Material</DataItem>
  <DataItem Name="Mass" Format="HDF" Datatype="Int"  Dimensions="%(nn)s">%(fgeom)s:/Mass</DataItem>
  <DataItem Name="Jac" Format="HDF" Datatype="Int"  Dimensions="%(nn)s">%(fgeom)s:/Jac</DataItem>
  <Grid Name="mesh.0.0">
    <Time Value="0.0"/>
    <Topology Type="Quadrilateral" NumberOfElements="%(ne)s">
    <DataItem Format="HDF" Datatype="Int" Dimensions="%(ne)s 4">
    %(fgeom)s:/Elements
    </DataItem>
    </Topology>
    <Geometry Type="XY">
    <DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="%(nn)s 2">
    %(fgeom)s:/Nodes
    </DataItem>
    </Geometry>
    <Attribute Name="Displ" Center="Node" AttributeType="Vector" Dimensions="%(nn)s 2">
    <DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="%(nn)s 2">
    %(fdata)s:/displ
    </DataItem>
    </Attribute>
    <Attribute Name="Veloc" Center="Node" AttributeType="Vector" Dimensions="%(nn)s 2">
    <DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="%(nn)s 2">
    %(fdata)s:/veloc
    </DataItem>
    </Attribute>
    <Attribute Name="Domain" Center="Grid" AttributeType="Scalar" Dimensions="1">
    <DataItem Format="XML" Datatype="Int"  Dimensions="1">0</DataItem>
    </Attribute>
    <Attribute Name="Mat" Center="Cell" AttributeType="Scalar" Dimensions="%(ne)s">
    <DataItem Reference="XML">/Xdmf/Domain/Grid/DataItem[@Name="Mat"]</DataItem>
    </Attribute>
    <Attribute Name="Mass" Center="Node" AttributeType="Scalar" Dimensions="%(nn)s">
    <DataItem Reference="XML">/Xdmf/Domain/Grid/DataItem[@Name="Mass"]</DataItem>
    </Attribute>
    <Attribute Name="Jac" Center="Node" AttributeType="Scalar" Dimensions="%(nn)s">
    <DataItem Reference="XML">/Xdmf/Domain/Grid/DataItem[@Name="Jac"]</DataItem>
    </Attribute>
  </Grid>
</Grid>
</Domain>
</Xdmf>
"""


def write_mesh_h5(fname, mesh, proc):
    NE = mesh.nelements()
    NN = mesh.nnodes()
    f = h5py.File(fname+".h5","w")
    f["Elements"] = mesh.elements
    f["Nodes"] = mesh.nodes
    f["displ"] = np.zeros((NN,2), float)
    f["veloc"] = np.zeros((NN,2), float)
    f["Material"] = np.zeros((NE,), int)
    f["Mass"] = np.zeros((NN,), int)
    f["Jac"] = np.zeros((NN,), int)
    g = file(fname+".xmf","w")
    d = {"fgeom" : fname+".h5",
         "fdata" : fname+".h5",
         "ne" : NE,
         "nn" : NN,
     }
    g.write(TMPL_H5 % d)



def gen_rect_mesh(x, y):
    X, Y = np.meshgrid(x,y)
    NX, NY = X.shape
    print NX, NY
    NN = NX*NY
    NE = (NX-1)*(NY-1)
    nodes = np.zeros( (NN,2), float)
    nodes[:,0] = X.flat
    nodes[:,1] = Y.flat

    mesh = Mesh2D(nodes)

    nid = np.arange(NN)
    nid.shape = (NX, NY)
    # elements to node mapping
    elements = np.zeros( (NX-1, NY-1, 4), int)
    elements[:,:,0] = nid[:-1,1:]
    elements[:,:,1] = nid[1:,1:]
    elements[:,:,2] = nid[1:,:-1]
    elements[:,:,3] = nid[:-1,:-1]
    elements.shape = -1, 4

    mesh.compute_connectivity(elements)
    # Material
    mesh.materials = np.zeros( (NE,), int)

    return mesh



# demo :
NX = 200
NY = 100
x = np.linspace(-3000, 3000, NX)
y = np.linspace(-3000,    0, NY)

mesh = gen_rect_mesh(x,y)
print mesh.nodes
print mesh.elements

write_mesh("mesh4spec.0000", mesh, 0)
write_mesh_h5("mesh2D", mesh, 0)
