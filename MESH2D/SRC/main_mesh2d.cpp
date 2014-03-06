
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#define OMPI_SKIP_MPICXX
#include <hdf5.h>
#include <cassert>
#include <unistd.h>
#include "metis.h"

using std::string;
using std::vector;
using std::pair;
using std::map;

struct Point {
    Point() {}
    Point(double u, double v):x(u),y(v) {}
    Point(const Point& p):x(p.x),y(p.y) {}
    double x,y;
};

struct Quad {
    int n[4];
};

struct Edge {
    int e[2];
};

// Edge index : <node1,node2> with node1<node2
typedef pair<int,int> edge_idx_t;

edge_idx_t make_edge_idx(int n1, int n2)
{
    if (n2<n1) return edge_idx_t(n2,n1);
    else return edge_idx_t(n1,n2);
}

struct edge_info_t
{
    edge_info_t():edge(-1),coherency(-1) {}
    edge_info_t(int e, int c):edge(e),coherency(c) {}
    edge_info_t(const edge_info_t& ei):edge(ei.edge),coherency(ei.coherency) {}
    int edge;
    int coherency;
};

struct Comm_proc {
    vector<int> m_vertices;
    map<int,int> m_vertices_map;
    vector<int> m_edges;
    map<edge_idx_t,edge_info_t> m_edges_map;
    vector<int> m_coherency;
};

struct MeshProcInfo {
    // Computed for storage
    vector<double> m_nodes;
    vector<int> m_quads; // 4 vertice number for each quad
    vector<int> m_material;
    vector<int> m_edges; // 4 edge number for each quad
    vector<int> m_edges_elems; // Elem0,Elem1 for each edge
    vector<int> m_edges_local; // local face number of the corresponding element (Which_face)
    vector<int> m_edges_vertices; // vert0,vert1 for each edge

    int n_elements() const { return m_quads.size()/4; }
    int n_edges() const { return m_edges_vertices.size()/2; }
    int n_vertices() const { return m_nodes.size()/2; }

    // Bookkeeping info
    vector<int> m_node_map;
    vector<int> m_quad_map;
    map<edge_idx_t,int> m_edge_map; // maps edge_idx to edge num+1 (so that 0 is used as not assigned)


    // Methods
    void add_edge_element(Quad& q, int fn, int en, int qn);
    void create_new_edge();
    void get_quad_edge(Quad& q, int fn, int& v0, int& v1);

    map<int,Comm_proc> m_comm;
};

class Mesh2D {
public:

    void read_mesh(const string& fname);
    void partition_metis(int nproc);
    void partition_scotch(int nproc);
    void check_cell_orient();
    void write_proc_field(const string& fname);
    void write_proc_file(const string& fname, int rk);
    void gather_proc_info(MeshProcInfo& info, int rk);
    void store_local_quad_points(MeshProcInfo& info, int qn, Quad& quad);
    void prepare_new_proc_info(MeshProcInfo& info);
    int m_nprocs;
    int m_mat_max;
    vector<double> m_px;
    vector<double> m_py;
    vector<Quad> m_quads;
    vector<int>  m_procs;
    vector<int>  m_mat1;
    vector<int>  m_mat2;
};


void Mesh2D::store_local_quad_points(MeshProcInfo& info, int qn, Quad& quad)
{
    // Create new local quad number
    info.m_quad_map[qn] = info.n_elements();
    for(int k=0;k<4;++k) {
	int nn = quad.n[k];
	int local_num = info.m_node_map[nn];
	if (local_num==-1) {
	    local_num = info.m_nodes.size()/2;
	    info.m_nodes.push_back(m_px[nn]);
	    info.m_nodes.push_back(m_py[nn]);
	    info.m_node_map[nn] = local_num;
	}
    }
}

void Mesh2D::prepare_new_proc_info(MeshProcInfo& info)
{
    info.m_node_map.clear();
    info.m_material.clear();
    info.m_nodes.clear();
    info.m_quad_map.clear();
    info.m_node_map.resize(m_px.size(), -1);
    info.m_quad_map.resize(m_quads.size(), -1);
}

void MeshProcInfo::get_quad_edge(Quad& q, int fn, int& v0, int& v1)
{
    switch(fn) {
    case 0:
	v0 = q.n[0];
	v1 = q.n[1];
	break;
    case 1:
	v0 = q.n[1];
	v1 = q.n[2];
	break;
    case 2:
	v0 = q.n[3];
	v1 = q.n[2];
	break;
    case 3:
	v0 = q.n[0];
	v1 = q.n[3];
	break;
    default:
	printf("ERR: internal error in get_quad_edge\n");
	exit(1);
    }
    v0 = m_node_map[v0];
    v1 = m_node_map[v1];
}

void MeshProcInfo::create_new_edge()
{
    m_edges_elems.push_back(-1);
    m_edges_elems.push_back(-1);
    m_edges_local.push_back(-1);
    m_edges_local.push_back(-1);
    m_edges_vertices.push_back(-1);
    m_edges_vertices.push_back(-1);
}

void MeshProcInfo::add_edge_element(Quad& q, int fn, int en, int qn)
{
    // fn : local face number of current quad
    // qn : quad number
    // en : edge number
    int j=0;
    int v0=-1, v1=-1;
    int lqn = m_quad_map[qn];
    assert(m_edges_elems[2*en+1]==-1);
    if (m_edges_elems[2*en]==-1) {
	j=0;
	get_quad_edge(q, fn, v0, v1);
    } else if (m_edges_elems[2*en+0]>lqn) {
	j=0;
	// Move elem0 to index 1
	m_edges_elems[2*en+1] = m_edges_elems[2*en+0];
	m_edges_local[2*en+1] = m_edges_local[2*en+0];
	// Reset vertices for elem0
	get_quad_edge(q, fn, v0, v1);
    } else {
	j=1;
    }
    m_edges_elems[2*en+j] = lqn;
    m_edges_local[2*en+j] = fn;
    if (j==0) {
	m_edges_vertices[2*en+0] = v0;
	m_edges_vertices[2*en+1] = v1;
    } else {
	v0 = m_edges_vertices[2*en+0];
	v1 = m_edges_vertices[2*en+1];
    }
    //printf("qn=%d fn=%d en=%d : j=%d v0=%d v1=%d\n", qn, fn, en, j, v0, v1);
}

void Mesh2D::gather_proc_info(MeshProcInfo& info, int rk)
{
    prepare_new_proc_info(info);
    for(int qn=0;qn<m_quads.size();++qn) {
	Quad& quad = m_quads[qn];
	if (rk!=m_procs[qn]) {
	    continue;
	}
	store_local_quad_points(info, qn, quad);
	info.m_material.push_back(m_mat1[qn]);
	info.m_material.push_back(0); // Flag solid/fluid... TODO remove
	info.m_material.push_back(m_mat2[qn]);
	for(int k=0;k<4;++k) {
	    info.m_quads.push_back(info.m_node_map[quad.n[k]]);
	    edge_idx_t e = make_edge_idx(quad.n[k], quad.n[(k+1)%4]);
	    int& edge_num = info.m_edge_map[e]; //edge_num is offset by one so that zero means unassigned
	    if (edge_num==0) {
		edge_num = info.m_edge_map.size();
		info.create_new_edge();
		assert(edge_num==info.m_edges_local.size()/2);
	    }
	    info.add_edge_element(quad, k, edge_num-1, qn);
	    info.m_edges.push_back(edge_num-1);
	}
    }
    for(int k=0;k<info.m_edges_local.size();k+=2) assert(info.m_edges_local[k]!=-1);
    for(int qn=0;qn<m_quads.size();++qn) {
	Quad& quad = m_quads[qn];
	if (rk!=m_procs[qn]) {
	    // mark edges and vertices involved in communications
	    for(int k=0;k<4;++k) {
		// Check vertices
		int local_node_num = info.m_node_map[quad.n[k]];
		if (local_node_num>=0) {
		    Comm_proc& comm = info.m_comm[m_procs[qn]];
		    comm.m_vertices_map[quad.n[k]] = local_node_num;
		}
		// Check edges
		edge_idx_t e = make_edge_idx(quad.n[k], quad.n[(k+1)%4]);
		int edge_num = info.m_edge_map[e]-1; //edge_num is offset by one so that zero means unassigned
		if (edge_num>=0) {
		    Comm_proc& comm = info.m_comm[m_procs[qn]];
		    comm.m_edges_map[e]=edge_info_t(edge_num,1);
		}
	    }
	}
    }
    map<int,Comm_proc>::iterator cit;
    for(cit=info.m_comm.begin();cit!=info.m_comm.end();++cit) {
	Comm_proc& comm=cit->second;
	map<edge_idx_t,edge_info_t>::const_iterator eit;
	for(eit=comm.m_edges_map.begin();eit!=comm.m_edges_map.end();++eit) {
	    comm.m_edges.push_back(eit->second.edge);
	    comm.m_coherency.push_back(eit->second.coherency);
	}
	map<int,int>::const_iterator vit;
	for(vit=comm.m_vertices_map.begin();vit!=comm.m_vertices_map.end();++vit) {
	    comm.m_vertices.push_back(vit->second);
	}
    }
}

int read_attr_int(hid_t dset_id, const char* attrname)
{
    hsize_t dims[1] = {1,};
    int res;
    hid_t attid = H5Aopen(dset_id, attrname, H5P_DEFAULT);
    //hid_t spcid = H5Aget_space(attid);
    //TODO check dims of attr
    H5Aread(attid, H5T_NATIVE_INT, &res);
    H5Aclose(attid);
    return res;
}

void write_attr_int(hid_t dset_id, const char* attrname, int val)
{
    hsize_t dims[1] = {1,};
    hid_t spcid = H5Screate(H5S_SCALAR);
    hid_t attid = H5Acreate2(dset_id, attrname, H5T_STD_I64LE, spcid, H5P_DEFAULT, H5P_DEFAULT);

    //TODO check dims of attr
    H5Awrite(attid, H5T_NATIVE_INT, &val);
    H5Aclose(attid);
    H5Sclose(spcid);
}

int get_dset1d_size(hid_t dset_id)
{
    hsize_t dims[1], maxdims[1];
    hid_t space_id = H5Dget_space(dset_id);
    int ndims = H5Sget_simple_extent_ndims(space_id);
    assert(ndims==1);
    H5Sget_simple_extent_dims(space_id, dims, maxdims);
    H5Sclose(space_id);
    return dims[0];
}

void get_dset2d_size(hid_t dset_id, int& d1, int& d2)
{
    hsize_t dims[2], maxdims[2];
    hid_t space_id = H5Dget_space(dset_id);
    int ndims = H5Sget_simple_extent_ndims(space_id);
    assert(ndims==2);
    H5Sget_simple_extent_dims(space_id, dims, maxdims);
    H5Sclose(space_id);
    d1 = dims[0];
    d2 = dims[1];
}

void read_dset_1d_d(hid_t g, const char* dname, vector<double>& v)
{
    hid_t dset_id;
    int dim;
    dset_id = H5Dopen2(g, dname, H5P_DEFAULT);
    dim = get_dset1d_size(dset_id);
    v.resize(dim);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
    H5Dclose(dset_id);
}

void read_dset_1d_i(hid_t g, const char* dname, vector<int>& v)
{
    hid_t dset_id;
    int dim;
    dset_id = H5Dopen2(g, dname, H5P_DEFAULT);
    dim = get_dset1d_size(dset_id);
    v.resize(dim);
    H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
    H5Dclose(dset_id);
}

void write_dset_1d_i(hid_t file_id, const char* dname, const vector<int>& v)
{
   hsize_t dims[1];
   herr_t  status;
   hid_t   dset_id;
   hid_t   dspc_id;

   if (H5Lexists(file_id, dname, H5P_DEFAULT)) {
       dset_id = H5Dopen(file_id, dname, H5P_DEFAULT);
       // TODO Check size coherency
    } else {
	/* Create the data space for the dataset. */
	dims[0] = v.size();
	dspc_id = H5Screate_simple(1, dims, NULL);
	dset_id = H5Dcreate2(file_id, dname, H5T_STD_I32LE, dspc_id,
				      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Sclose(dspc_id);
   }
   H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
   status = H5Dclose(dset_id);
}

void write_dset_2d_i(hid_t file_id, const char* dname, int d2, const vector<int>& v)
{
   hsize_t     dims[2];
   herr_t      status;

   /* Create the data space for the dataset. */
   dims[0] = v.size()/d2;
   dims[1] = d2;
   hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
   hid_t dataset_id = H5Dcreate2(file_id, dname, H5T_STD_I32LE, dataspace_id,
				 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
   status = H5Dclose(dataset_id);
   status = H5Sclose(dataspace_id);
}

void write_dset_2d_r(hid_t file_id, const char* dname, int d2, const vector<double>& v)
{
   hsize_t     dims[2];
   herr_t      status;

   /* Create the data space for the dataset. */
   dims[0] = v.size()/d2;
   dims[1] = d2;
   hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
   hid_t dataset_id = H5Dcreate2(file_id, dname, H5T_IEEE_F64LE, dataspace_id,
				 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
   status = H5Dclose(dataset_id);
   status = H5Sclose(dataspace_id);
}

void read_quads(hid_t g, const char* dname, vector<Quad>& v)
{
    hid_t dset_id;
    int n0, n1;
    dset_id = H5Dopen2(g, dname, H5P_DEFAULT);
    get_dset2d_size(dset_id, n0, n1);

    assert(n1==4);
    v.resize(n0);
    H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
    H5Dclose(dset_id);
}

void Mesh2D::read_mesh(const string& fname)
{
    hid_t file_id;

    file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    // Read points
    read_dset_1d_d(file_id, "X", m_px);
    read_dset_1d_d(file_id, "Y", m_py);
    read_quads(file_id, "Q", m_quads);
    read_dset_1d_i(file_id, "M1", m_mat1);
    read_dset_1d_i(file_id, "M2", m_mat2);
    assert(m_mat1.size()==m_mat2.size());
    assert(m_quads.size()==m_mat1.size());
    printf("%ld nodes, %ld Quads\n", m_px.size(), m_quads.size());
    H5Fclose(file_id);
    m_nprocs = 1;
    m_procs.resize(m_quads.size(), 0);
    m_mat_max = 0;
    for(int k=0;k<m_mat1.size();++k) if (m_mat1[k]>m_mat_max) m_mat_max = m_mat1[k];

}

void Mesh2D::check_cell_orient()
{
}

void Mesh2D::partition_metis(int nproc)
{
    idx_t options[METIS_NOPTIONS];
    idx_t ne = m_quads.size();
    idx_t nn = m_px.size();
    idx_t ncommon=1;
    idx_t numflag=0;
    vector<idx_t> eptr, eind;
    idx_t* xadj;
    idx_t* adjncy;

    m_nprocs = nproc;
    METIS_SetDefaultOptions(options);
    for(int k=0;k<METIS_NOPTIONS;++k) printf("OPT:%4d = %d\n", k, options[k]);
    options[METIS_OPTION_NUMBERING] = 0;
    //options[METIS_OPTION_CONTIG] = 1;
    //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

    // Setup graph for metis
    eptr.resize(ne+1);
    eind.resize(4*ne);
    for(int i=0;i<ne;++i) {
	eptr[i] = 4*i;
	for(int k=0;k<4;++k) eind[4*i+k] = m_quads[i].n[k];
    }
    eptr[ne]=4*ne;
    METIS_MeshToDual(&ne, &nn, &eptr[0], &eind[0], &ncommon, &numflag, &xadj, &adjncy);
    m_procs.resize(ne, -1);
    vector<idx_t> vwgt(ne, 1);
    vector<idx_t> vsize(ne, 1);
    vector<idx_t> adjwgt(xadj[ne],1);
    vector<real_t> tpwgts(nproc);
    idx_t edgecut;
    idx_t ncon=1;
    idx_t inproc = nproc;
    real_t ubvec[1] = { 1.001 }; // load imbalance
    for(int k=0;k<nproc;++k) tpwgts[k] = 1./nproc;
    METIS_PartGraphRecursive(&ne, &ncon, xadj, adjncy,
			     &vwgt[0], &vsize[0], &adjwgt[0], &inproc,
			     &tpwgts[0], ubvec, options, &edgecut, &m_procs[0]);
//    for(int k=0;k<ne;++k) printf("%d : %d\n", k, m_procs[k]);
}

void Mesh2D::write_proc_field(const string& fname)
{
    hid_t file_id;

    file_id = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Read points
    write_dset_1d_i(file_id, "Proc", m_procs);
    H5Fclose(file_id);
}


void Mesh2D::write_proc_file(const string& fname, int rk)
{
    hid_t fid;
    if (access(fname.c_str(), F_OK)==0) {
	// File already exists
	unlink(fname.c_str());
    }
    fid = H5Fopen(fname.c_str(), H5F_ACC_RDWR|H5F_ACC_CREAT, H5P_DEFAULT);

    MeshProcInfo info;

    gather_proc_info(info, rk);
    write_attr_int(fid, "ndim", 2);
    //m_nprocs=1;
    write_attr_int(fid, "n_processors", m_nprocs);
    write_attr_int(fid, "n_materials", m_mat_max+1);
    write_attr_int(fid, "n_elements", info.n_elements());
    write_attr_int(fid, "n_edges", info.n_edges());
    write_attr_int(fid, "n_vertices", info.n_vertices());
    write_dset_2d_r(fid, "nodes", 2, info.m_nodes);
    write_dset_2d_i(fid, "material", 3, info.m_material);
    write_dset_2d_i(fid, "elements", 4, info.m_quads);
    write_dset_2d_i(fid, "edges", 4, info.m_edges);
    // With linear Quad, vertices==elements
    write_dset_2d_i(fid, "vertices", 4, info.m_quads);
    write_dset_2d_i(fid, "faces_elem", 2, info.m_edges_elems);
    write_dset_2d_i(fid, "faces_which", 2, info.m_edges_local);
    write_dset_2d_i(fid, "faces_vertex", 2, info.m_edges_vertices);

    int n_comm = info.m_comm.size();
    //n_comm = 0;
    write_attr_int(fid, "n_communications", n_comm);
    map<int,Comm_proc>::const_iterator it;
    int comm_count=0;
    for(it=info.m_comm.begin();it!=info.m_comm.end();++it) {
	char grp_name[60];
	snprintf(grp_name, 60, "Comm%05d", comm_count);
	const Comm_proc& comm = it->second;
	hid_t grp = H5Gcreate(fid, grp_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	write_attr_int(grp, "processor", it->first);
	write_dset_1d_i(grp, "vertices", comm.m_vertices);
	write_dset_1d_i(grp, "edges", comm.m_edges);
	write_dset_1d_i(grp, "coherency", comm.m_coherency);
	H5Gclose(grp);
	comm_count++;
    }

    H5Fclose(fid);
}

int main(int argc, char** argv)
{
    char fname[1024];
    if (argc<4) {
	printf("Usage: mesh2dc Nproc mesh_input.h5 base_out\n");
	exit(1);
    }

    int nproc = atoi(argv[1]);
    string fmesh = argv[2];

    Mesh2D mesh;

    mesh.read_mesh(fmesh);
    if (nproc>1) {
	mesh.partition_metis(nproc);
    } else {
    }
    mesh.write_proc_field(fmesh);

    for(int k=0;k<nproc;++k) {
	snprintf(fname, 1024,"%s.%04d.h5", argv[3], k);
	mesh.write_proc_file(fname, k);
    }
}
// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=4 et tw=80 smartindent :*/
