/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include <cstdio>
#include "mesh.h"
#include "metis.h"
#include <map>
#include <cstdlib>

using std::map;
using std::multimap;
using std::vector;
using std::pair;

// =====================================================


/// Reference numbering for faces
static FaceDesc RefFace[6] = {
    { {0,1,2,3}, {0, 1, 2, 3} },
    { {0,1,5,4}, {0, 4, 5, 6} },
    { {1,2,6,5}, {1, 7, 8, 4} },
    { {3,2,6,7}, {2, 7, 9,10} },
    { {0,3,7,4}, {3,10,11, 6} },
    { {4,5,6,7}, {5, 8, 9,11} },
};

/// Reference numbering of edges
static int RefEdge[12][2] = {
    { 0, 1},
    { 1, 2},
    { 3, 2},
    { 0, 3},
    { 1, 5},
    { 4, 5},
    { 0, 4},
    { 2, 6},
    { 5, 6},
    { 7, 6},
    { 3, 7},
    { 4, 7}
};

/// List the faces each node belongs to
static int NodeToFace[8][3] = {
    {0,1,4}, // 0
    {0,1,2}, // 1
    {0,2,3}, // 2
    {0,3,4}, // 3
    {1,4,5}, // 4
    {1,2,5}, // 5
    {2,3,5}, // 6
    {3,4,5}, // 7
};


// Given orientation provide the permutation that reorders the nodes of the permuted
// face into face of orientation 0
static int Face_Orient[8][4] = {
    {0,1,2,3},
    {1,0,3,2},
    {3,2,1,0},
    {2,3,0,1},
    {0,3,2,1},
    {3,0,1,2},
    {1,2,3,0},
    {2,1,0,3},
};

// Returns an int representing face orientation relative to reference face
// 0: 0,1,2,3 : reference
// 6: 1,2,3,0
// 3: 2,3,0,1
// 5: 3,0,1,2
// 4: 0,3,2,1
// 2: 3,2,1,0
// 7: 2,1,0,3
// 1: 1,0,3,2
static int face_orientation( int nf, FaceDesc& f )
{
    FaceDesc& r = RefFace[nf];

    if (f.v[0] == r.v[0] && f.v[1] == r.v[1] && f.v[2] == r.v[2] && f.v[3] == r.v[3]) return 0;
    if (f.v[0] == r.v[1] && f.v[1] == r.v[2] && f.v[2] == r.v[3] && f.v[3] == r.v[0]) return 6;
    if (f.v[0] == r.v[2] && f.v[1] == r.v[3] && f.v[2] == r.v[0] && f.v[3] == r.v[1]) return 3;
    if (f.v[0] == r.v[3] && f.v[1] == r.v[0] && f.v[2] == r.v[1] && f.v[3] == r.v[2]) return 5;
    if (f.v[0] == r.v[0] && f.v[1] == r.v[3] && f.v[2] == r.v[2] && f.v[3] == r.v[1]) return 4;
    if (f.v[0] == r.v[3] && f.v[1] == r.v[2] && f.v[2] == r.v[1] && f.v[3] == r.v[0]) return 2;
    if (f.v[0] == r.v[2] && f.v[1] == r.v[1] && f.v[2] == r.v[0] && f.v[3] == r.v[3]) return 7;
    if (f.v[0] == r.v[1] && f.v[1] == r.v[0] && f.v[2] == r.v[3] && f.v[3] == r.v[2]) return 1;
    // The face is probably twisted...
    return -1;
}

static int find_face(FaceDesc& f)
{
    int fc[3];
    int n,v,k,l,u;
    int c=0;
    n = f.v[0];
    // We setup the list of the 3 possible face that vertex 0 belongs to
    for(k=0;k<3;++k) fc[k] = NodeToFace[n][k];
    // And we proceed to eliminate faces that mismatch (ie those from vertex 0
    // that do not belong to the list of possible faces of other vertices)
    for(int v=1;v<4;++v) {
	//printf("a.%d.%d.%d\n",fc[0],fc[1],fc[2]);
	n = f.v[v];
	for(k=0;k<3;++k) {
	    if (fc[k]==-1) continue;
	    u = k;
	    for(l=0;l<3;++l) {
		if (fc[k]==NodeToFace[n][l]) break;
	    }
	    if (l==3) {
		fc[k]=-1;
		//printf("b.%d.%d.%d\n",fc[0],fc[1],fc[2]);
		c++;
		// Reaching here means we eliminated all faces. Bad.
		if (c==3) return -1;
	    }
	}
    }
    // There should be only one face left
    if (c!=2) return -1;
    // And it should be fc[u], u being the last index visited in the above loop
    return fc[u];
}

static int which_face(FaceDesc& fc, int orient)
{
    int nf,k;
    int oriented_face[4];
    for(k=0;k<4;++k)
	oriented_face[k] = fc.v[Face_Orient[orient][k]];
    for(nf=0;nf<6;++nf) {
	for(k=0;k<4;++k) {
	    if (oriented_face[k]!=RefFace[nf].v[k]) break;
	}
	if (k==4) break;
    }
    return nf;
}

bool Mesh3D::shared_face(int el0, int nf0, int el1, FaceDesc& other)
{
    const FaceDesc& fc = RefFace[nf0];
    int otherface[4] = {-1,-1,-1,-1};
    int k;
    int i0 = m_elems_offs[el0];
    int i1 = m_elems_offs[el1];
    for(int n=0;n<4;++n) {
	int nn = m_elems[i0+fc.v[n]];
	for(k=0;k<8;++k) {
	    if (nn == m_elems[i1+k]) {
		other.v[n] = k;
		break;
	    }
	}
	if (k==8) return false;
    }
    return true;
}



// =====================================================
int Mesh3D::add_node(double x, double y, double z)
{
    m_xco.push_back( x );
    m_yco.push_back( y );
    m_zco.push_back( z );
    return m_xco.size()-1;
}

int Mesh3D::add_elem(int mat_idx, const HexElem& el)
{
    // Builds elem<->vertex graph
    for(int i=0;i<8;++i) {
	m_elems.push_back(el.v[i]);
    }
    m_elems_offs.push_back(m_elems.size());
    m_mat.push_back( mat_idx );
}


void Mesh3D::partition_mesh(int n_parts)
{
    int ne = n_elems();
    int nn = n_vertices();
    int ncommon = 1;
    int numflags = 0;
    int ncon=1;
    int *vwgt=0L;
    int *vsize=0L;
    int *adjwgt=0L;
    float *tpwgts=0L;
    float *ubvec=0L;
    int *options=0L;
    int edgecut;

    n_procs = n_parts;
    m_procs.resize(ne);

    METIS_MeshToDual(&ne, &nn, &m_elems_offs[0], &m_elems[0],
		     &ncommon, &numflags, &m_xadj, &m_adjncy);

    METIS_PartGraphKway(&ne, &ncon, m_xadj, m_adjncy,
			vwgt, vsize, adjwgt, &n_procs, tpwgts, ubvec,
			options, &edgecut, &m_procs[0]);

    compute_comm_elements();
}

void Mesh3D::compute_comm_elements()
{
    // Note for future parallelisation : in this routine, each part
    // only needs to know about it's direct neighbours' elements
    // So that a parallel algorithm could provide a subgraph to each
    // processor containing elements of a part and the ghost elements.

    // Find adjacent elements that resides on different processors
    // Note: the graph will be symetric.
    for(int gel=0;gel<n_elems();++gel) {
	for(int i=m_xadj[gel];i<m_xadj[gel+1];++i) {
	    int nel = m_adjncy[i];
	    if (m_procs[gel]!=m_procs[nel]) {
		pair<int,int> proc_pair(gel, nel);
		m_elem_proc_map.insert( proc_pair );
		// XXX: avoid copy if already exists
		m_shared[proc_pair] = shared_part_t();
	    }
	}
    }
    // The shared info has to be exchanged at some point so that
    // its coherent . There should be a synchronisation point where
    // proc p1>p0 waits for p0 to send it's shared_part_t_ data
}

void Mesh3D::compute_shared_objects_for_part(int part)
{
    multimap<int,int>::const_iterator it;
    // Finds faces, edges, vertices shared, and their relative orientation
    // the element on the lowest rank processor receives orientation 0
    for(it=m_elem_proc_map.begin();it!=m_elem_proc_map.end();++it) {
	int local_el, other_el;
	local_el = it->first;
	if (m_procs[local_el]!=part) continue;
	other_el = it->second;
	// Faces : we already have them from add_local_face
    }
}

void Mesh3D::compute_local_part(int part, MeshPart& loc)
{
    loc.part = part;

    loc.n_elems_per_proc = 0;
    loc.m_n_faces = 0;
    loc.m_n_edges = 0;
    // Count number of elements on this proc
    loc.e_g2l.resize(n_elems(),-1);
    for(int i=0;i<n_elems();++i) {
	if (elem_part(i) == loc.part) {
	    loc.e_l2g.push_back(i);
	    loc.e_g2l[i] = loc.n_elems_per_proc;
	    loc.n_elems_per_proc++;
	}
    }
    compute_nodes_indexes(loc);
    compute_local_connectivity(loc);
}

void Mesh3D::compute_nodes_indexes(MeshPart& loc)
{
    loc.n_g2l.resize(n_nodes(), -1);

    for(int i=0;i<n_elems();++i) {
	if (elem_part(i)!=loc.part) continue;
	for(int k=m_elems_offs[i];k<m_elems_offs[i+1];++k) {
	    int node = m_elems[k];
	    if (loc.n_g2l[node]==-1) {
		// First visit, assign node local number
		loc.n_g2l[node] = loc.n_l2g.size();
		loc.n_l2g.push_back(node);
	    }
	}
    }
}

void Mesh3D::add_local_face(MeshPart& loc, int iel, int gel, int nf)
{
    const FaceDesc& fc = RefFace[nf];
    // For each neighbour :
    int found=0;
    FaceDesc other;
    int neighbor;
    for(int i=m_xadj[gel];i<m_xadj[gel+1];++i) {
	neighbor = m_adjncy[i];
	if (shared_face(gel, nf, neighbor, other)) {
	    found = 1;
	    break;
	}
    }
    if (found==0) {
	// The face is not shared
	loc.faces[6*iel+nf] = loc.m_n_faces;
	loc.faces_orient[6*iel+nf] = 0;
	loc.m_n_faces++;
	return;
    }

    // Found a neighbor sharing a face:
    if (elem_part(neighbor) == loc.part) {
	// Both elements are on the same processor
	// Since local numbers are assigned following increasing
	// global elements indexes :
	if (neighbor>gel) {
	    // Neighbor hasn't been visited yet
	    loc.faces[6*iel+nf] = loc.m_n_faces;
	    loc.faces_orient[6*iel+nf] = 0;
	    //printf("New:%d o=0 n=%d\n", nf, loc.m_n_faces);
	    loc.m_n_faces++;
	} else {
	    //
	    int onf = find_face(other);
	    int orientation = face_orientation(onf, other);
	    int lk = loc.e_g2l[neighbor];
	    loc.faces[6*iel+nf] = loc.faces[6*lk+onf];
	    loc.faces_orient[6*iel+nf] = orientation;
	    //other.show_face();
	    //printf("Face:%d o=%d n=%d\n", onf, orientation, loc.faces[6*iel+nf]);
	}
    } else {
	// Not on the same processor
	loc.faces[6*iel+nf] = loc.m_n_faces;
	loc.faces_orient[6*iel+nf] = 0;
	int num = elem_part(neighbor);
	// TODO: register previous proc face number
	if (num<loc.part) {
	    int onf = find_face(other);
	    int orientation = face_orientation(onf, other);
	    //m_shared[num]
	} else {
	}
	loc.m_n_faces++;
    }
}

void Mesh3D::add_local_edge(MeshPart& loc, int iel, int gel, int ne)
{
    int ve0, ve1; // vertex nums of edge

    ve0 = m_elems[8*gel+RefEdge[ne][0]];
    ve1 = m_elems[8*gel+RefEdge[ne][1]];
    ordered_edge_t edg(ve0,ve1);

    map<ordered_edge_t,edge_t>::iterator it;

    it = m_edge_map.find(edg);

    if (it==m_edge_map.end()) {
	// Edge never visited
	// New number
	m_edge_map[edg] = edge_t(ve0, ve1, loc.m_n_edges, gel);
	loc.edges[12*iel+ne] = loc.m_n_edges;
	loc.edges_orient[12*iel+ne] = 0;
	loc.m_n_edges++;
    } else {
	edge_t& edg0 = it->second;
	loc.edges[12*iel+ne] = edg0.n;
	loc.edges_orient[12*iel+ne] = edg0.orient(ve0,ve1);
	edg0.elems.push_back(gel);
    }
}

void Mesh3D::compute_local_connectivity(MeshPart& loc)
{
    loc.m_n_faces = 0;

    loc.faces.resize(6*loc.n_elems(),-1);
    loc.faces_orient.resize(6*loc.n_elems_per_proc,-1);
    loc.edges.resize(12*loc.n_elems(),-1);
    loc.edges_orient.resize(12*loc.n_elems(),-1);
    loc.vertices.resize(8*loc.n_elems(),-1);
    for(int iel=0;iel<loc.n_elems_per_proc;++iel) {
	int gel = loc.e_l2g[iel];
	for(int nf=0;nf<6;++nf) {
	    add_local_face(loc, iel, gel, nf);
	}
	for(int ne=0;ne<11;++ne) {
	    add_local_edge(loc, iel, gel, ne);
	}
    }
}

int Mesh3D::read_materials(const std::string& str)
{
    FILE* f = fopen(str.c_str(), "r");
    int nmats;
    char type;
    double vs, vp, rho;
    int ngllx, nglly, ngllz;
    double dt, Qp, Qmu;

    fscanf(f, "%d", &nmats);
    for(int k=0;k<nmats;++k) {
        fscanf(f, "%c %lf %lf %lf %d %d %d %lf %lf %lf",
               &type, &vs, &vp, &rho, &ngllx, &nglly, &ngllz,
               &dt, &Qp, &Qmu);
        m_materials.push_back(Material(type, vp, vs, rho, Qp, Qmu, ngllx, nglly, ngllz));
    }
}

int Mesh3D::add_material()
{
    m_materials.push_back(Material());
}

void Mesh3D::get_local_material(MeshPart& loc, std::vector<int>& imat)
{
    imat.clear();
    for(int k=0; k<loc.n_elems(); ++k) {
	int g = loc.e_l2g[k];
	imat.push_back(m_mat[g]);
	imat.push_back( (m_materials[m_mat[g]].is_fluid() ? 0 : 1) );
    }
}

void Mesh3D::get_local_nodes(MeshPart& loc, std::vector<double>& tmp)
{
    tmp.clear();
    for(int k=0; k<loc.n_nodes();++k) {
	int g = loc.n_l2g[k];
	tmp.push_back(m_xco[g]);
	tmp.push_back(m_yco[g]);
	tmp.push_back(m_zco[g]);
    }
}

void Mesh3D::get_local_faces(MeshPart& loc, std::vector<int>& tmp)
{
    tmp.clear();
}

void Mesh3D::get_local_elements(MeshPart& loc, std::vector<int>& tmp)
{
    tmp.clear();
    for(int k=0;k<loc.n_elems();++k) {
	int g = loc.e_l2g[k];
	for(int i=m_elems_offs[g];i<m_elems_offs[g+1];++i) {
	    int ln = loc.n_g2l[m_elems[i]];
	    tmp.push_back(ln);
	}
    }
}

void Mesh3D::read_mesh_file(const std::string& fname)
{
    hid_t file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    h5h_read_dset_Nx3(file_id, "/Nodes", m_xco, m_yco, m_zco);

    hid_t dset_id;
    hsize_t n0, n1;
    if (H5Lexists(file_id, "/Sem3D/Hexa8", H5P_DEFAULT)) {
        read_mesh_hexa8(file_id);
    } else if (H5Lexists(file_id, "/Sem3D/Hexa27", H5P_DEFAULT)) {
        read_mesh_hexa27(file_id);
    }
    else {
        printf("ERR: only Quad4 and Quad8 are supported\n");
        exit(1);
    }

    h5h_read_dset(file_id, "/Sem3D/Mat", m_mat);
}

void Mesh3D::read_mesh_hexa8(hid_t file_id)
{
    int nel, nnodes;
    h5h_read_dset_2d(file_id, "/Sem3D/Hexa8", nel, nnodes, m_elems);
    n_ctl_nodes = 8;
    if (nnodes!=8) {
        printf("Error: dataset /Sem/Hexa8 is not of size NEL*8\n");
        exit(1);
    }
    for(int k=0;k<nel;++k) {
        m_elems_offs.push_back(8*(k+1));
    }
}

void Mesh3D::read_mesh_hexa27(hid_t file_id)
{
}

void Mesh3D::build_vertex_to_elem_map()
{
    int nel = n_elems();
    m_vertex_to_elem.init(nel);
    for(int i=0;i<nel;++i) {
        for(int k=m_elems_offs[i];k<m_elems_offs[i+1];++k) {
            m_vertex_to_elem.add_link(m_elems[k], i);
        }
    }
}


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
