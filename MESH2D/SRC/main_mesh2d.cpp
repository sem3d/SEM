
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

class Point {
    public:
        Point() {}
        Point(double u, double v):x(u),y(v) {}
        Point(const Point& p):x(p.x),y(p.y) {}
        double x,y;
};

// Edge index : <node1,node2> with node1<node2
typedef pair<int,int> edge_idx_t;

class Quad {
    public:
        Quad() {nn = 0; n = NULL;}
        virtual ~Quad() {if(n){delete [] n; n = NULL; nn = 0;};}
        Quad(const Quad& Q) {nn = 0; n = NULL; *this = Q;}
        Quad& operator=(const Quad & Q) {
            if (this != &Q) {
                if(n)   {delete [] n; n = NULL; nn = 0;};
                if(Q.n) {nn = Q.nn; n = new int[nn]; for (int i=0; i<nn; i++) n[i] = Q.n[i];}
            }
            return *this;
        }
        int get_nb_nodes() {return nn;}
        int get_node_id(int i) {return (n && i>=0 && i<=nn-1) ? n[i] : -1;}
        void check_orient(vector<double>& x, vector<double>& y) {
            if (nn<=0 || !n) return;
            int n0=n[0]; assert(n0>=0 && n0<x.size());
            int n1=n[1]; assert(n1>=0 && n1<x.size());
            int n3=n[3]; assert(n3>=0 && n3<x.size());
            double v01[2] = {x[n1]-x[n0], y[n1]-y[n0]};
            double v03[2] = {x[n3]-x[n0], y[n3]-y[n0]};
            if ((v01[0]*v03[1]-v01[1]*v03[0])<0.){swap_orient(); printf("one quad (%f %f) has been swaped\n", x[n0]+v01[0]/2., y[n0]+v03[1]/2.);}
        };
        virtual void swap_orient() = 0;
        virtual edge_idx_t get_edge_from_node(int i) = 0; // Assume a node belongs to the upcoming edge (or current edge for intermediate node)
        virtual int get_face_from_node(int i) = 0; // Assume a node belongs to the upcoming face (or current edge for intermediate node)
        virtual bool is_intermediate_node(int i) = 0;
    protected:
        int nn;
        int* n;
};
class Quad4 : public Quad {
    public:
        Quad4() {nn = 0; n = NULL;}
        Quad4(const int* q):Quad() {nn = 4; n = new int[nn]; for(int i=0; q && i<nn; i++) n[i]=q[i];}
        virtual edge_idx_t get_edge_from_node(int i) // Return edge oriented the same way (low -> high) whatever element edge orientation may be
        {
            assert (i>=0 && i<=3);
            int ii=(i+1)%4; // In case i=3, ii=0
            return (n[i]<n[ii]) ? edge_idx_t(n[i],n[ii]) : edge_idx_t(n[ii],n[i]); // Reorder: low ID -> high ID
        };
        virtual void swap_orient() {int tmp=n[1]; n[1]=n[3]; n[3]=tmp;};
        virtual int get_face_from_node(int i) {assert (i>=0 && i<=3); return i;};
        virtual bool is_intermediate_node(int i) {assert (i>=0 && i<=3); return false;};
};
class Quad8 : public Quad4 { // Assume intermediate nodes are stored after principal nodes
    public:
        Quad8() {nn = 0; n = NULL;}
        Quad8(const int* q):Quad4() {nn = 8; n = new int[nn]; for(int i=0; q && i<nn; i++) n[i]=q[i];}
        virtual void swap_orient() {int tmp=n[4]; n[4]=n[7]; n[7]=tmp; tmp=n[5]; n[5]=n[6]; n[6]=tmp; Quad4::swap_orient();};
        virtual edge_idx_t get_edge_from_node(int i) {assert (i>=0 && i<=7); int ii = (i >= 4) ? i - 4 : i; return Quad4::get_edge_from_node(ii);};
        virtual int get_face_from_node(int i) {assert (i>=0 && i<=7); int ii = (i >= 4) ? i - 4 : i; return ii;};
        virtual bool is_intermediate_node(int i) {assert (i>=0 && i<=7); return (i >= 4) ? true : false;};
};

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

class MeshProcInfo {
    public:
        MeshProcInfo(const int npq) {m_npq = npq;}
        // Computed for storage
        vector<double> m_nodes;
        vector<int> m_quadnodes; // 4 or 8 node number for each quad
        vector<int> m_quadvertices; // always 4 vertice number for each quad
        vector<int> m_material;
        vector<int> m_edges; // 4 edge number for each quad
        vector<int> m_edges_elems; // Elem0,Elem1 for each edge
        vector<int> m_edges_local; // local face number of the corresponding element (Which_face)
        vector<int> m_edges_vertices; // vert0,vert1 for each edge

        int n_elements() const { return m_quadnodes.size()/m_npq; }
        int n_edges() const { return m_edges_vertices.size()/2; } // Divide by 2 as v0 and v1 are stored for each edge
        int n_vertices() const { return m_vert_map.size(); }

        // Bookkeeping info
        vector<int> m_node_map;
        map<int, int> m_vert_map; // for quad8, vertex global numbering != node numbering (intermediate nodes not considered)
        vector<int> m_quad_map;
        map<edge_idx_t,int> m_edge_map; // maps edge_idx to edge num+1 (so that 0 is used as not assigned)

        // Methods
        void add_edge_element(Quad& q, int fn, int en, int qn);
        void add_vertex(const int nid) {if(m_vert_map.find(nid)==m_vert_map.end()) {int vid=n_vertices(); m_vert_map[nid]=vid;}}
        void create_new_edge();
        void get_quad_edge(Quad& q, int fn, int& v0, int& v1);

        map<int,Comm_proc> m_comm;
        int m_npq; // Nb nodes per quad
};

class Mesh2D {
public:
    Mesh2D() {};
    virtual ~Mesh2D() {for(vector<Quad*>::iterator it = m_quads.begin(); it != m_quads.end(); ++it) {Quad* q=*it; if(q){delete q; q=NULL;}}};
    void read_mesh(const string& fname, int& t);
    void partition_metis(int nproc);
    void partition_scotch(int nproc);
    void check_cell_orient();
    void write_proc_field(const string& fname);
    void write_proc_file(const string& fname, const int t, int rk);
    void gather_proc_info(MeshProcInfo& info, int rk);
    void store_local_quad_points(MeshProcInfo& info, int qn, Quad& quad);
    void prepare_new_proc_info(MeshProcInfo& info);
    int m_nprocs;
    int m_mat_max;
    vector<double> m_px;
    vector<double> m_py;
    vector<Quad*> m_quads;
    vector<int>  m_procs;
    vector<int>  m_mat1;
    vector<int>  m_mat2;
};


void Mesh2D::store_local_quad_points(MeshProcInfo& info, int qn, Quad& quad)
{
    // Create new local quad number
    info.m_quad_map[qn] = info.n_elements();
    for(int k=0;k<quad.get_nb_nodes();++k) {
        int nn = quad.get_node_id(k);
        int local_num = info.m_node_map[nn];
        if (local_num==-1) {
            local_num = info.m_nodes.size()/2; // Divide by 2 as X and Y are stored for each node
            info.m_nodes.push_back(m_px[nn]);
            info.m_nodes.push_back(m_py[nn]);
            info.m_node_map[nn] = local_num;
            if(!quad.is_intermediate_node(k)) {
              info.add_vertex(local_num);
            }
        }
    }
}

void Mesh2D::prepare_new_proc_info(MeshProcInfo& info)
{
    info.m_node_map.clear();
    info.m_vert_map.clear();
    info.m_material.clear();
    info.m_nodes.clear();
    info.m_quad_map.clear();
    info.m_node_map.resize(m_px.size(), -1);
    info.m_quad_map.resize(m_quads.size(), -1);
}

void MeshProcInfo::get_quad_edge(Quad& q, int fn, int& v0, int& v1)
{
    int n0=-1, n1=-1;
    switch(fn) {
        case 0:
            n0 = q.get_node_id(0);
            n1 = q.get_node_id(1);
            break;
        case 1:
            n0 = q.get_node_id(1);
            n1 = q.get_node_id(2);
            break;
        case 2:
            n0 = q.get_node_id(3);
            n1 = q.get_node_id(2);
            break;
        case 3:
            n0 = q.get_node_id(0);
            n1 = q.get_node_id(3);
            break;
        default:
            printf("ERR: internal error in get_quad_edge\n");
            exit(1);
    }
    n0 = m_node_map[n0];
    n1 = m_node_map[n1];
    v0=m_vert_map[n0];
    v1=m_vert_map[n1];
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
    // Asserts we don't come here more than twice
    assert(m_edges_elems[2*en+1]==-1);
    if (m_edges_elems[2*en+0]==-1) {
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
        // Do nothing the vertices stored are already those of m_edges_elems[2*en+0]
        v0 = m_edges_vertices[2*en+0];
        v1 = m_edges_vertices[2*en+1];
    }
    //printf("qn=%d fn=%d en=%d : j=%d v0=%d v1=%d\n", qn, fn, en, j, v0, v1);
}

void Mesh2D::gather_proc_info(MeshProcInfo& info, int rk)
{
    prepare_new_proc_info(info);
    for(int qn=0;qn<m_quads.size();++qn) {
        Quad* quad = m_quads[qn];
        assert(quad);
        if (rk!=m_procs[qn]) {
            continue;
        }
        store_local_quad_points(info, qn, *quad);
        info.m_material.push_back(m_mat1[qn]);
        info.m_material.push_back(0); // flag solid/fluid... TODO remove
        info.m_material.push_back(m_mat2[qn]);
        for(int k=0;k<quad->get_nb_nodes();++k) {
            int locnode = info.m_node_map[quad->get_node_id(k)];
            info.m_quadnodes.push_back(locnode);
            if (quad->is_intermediate_node(k)) continue; // For topology construction, rely on principal nodes only
            info.m_quadvertices.push_back(info.m_vert_map[locnode]); // Log node as a vertex
            edge_idx_t e = quad->get_edge_from_node(k);
            int& edge_num = info.m_edge_map[e]; // edge_num is offset by one so that zero means unassigned
            if (edge_num==0) {
                edge_num = info.m_edge_map.size();
                info.create_new_edge();
                assert(edge_num==info.m_edges_local.size()/2);
            }
            info.add_edge_element(*quad, quad->get_face_from_node(k), edge_num-1, qn);
            info.m_edges.push_back(edge_num-1);
        }
    }
    for(int k=0;k<info.m_edges_local.size();k+=2) assert(info.m_edges_local[k]!=-1);
    for(int qn=0;qn<m_quads.size();++qn) {
        Quad* quad = m_quads[qn];
        assert(quad);
        if (rk!=m_procs[qn]) {
            // mark edges and vertices involved in communications
            for(int k=0;k<quad->get_nb_nodes();++k) {
                if (quad->is_intermediate_node(k)) continue; // For topology construction, rely on principal nodes only
                // Check vertices
                int local_node_num = info.m_node_map[quad->get_node_id(k)];
                if (local_node_num>=0) {
                    Comm_proc& comm = info.m_comm[m_procs[qn]];
                    comm.m_vertices_map[quad->get_node_id(k)] = local_node_num;
                }
                // Check edges
                edge_idx_t e = quad->get_edge_from_node(k);
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

void get_dset2d_size(hid_t dset_id, hsize_t& d1, hsize_t& d2)
{
    hsize_t dims[2], maxdims[2];
    hid_t space_id = H5Dget_space(dset_id);
    int ndims = H5Sget_simple_extent_ndims(space_id);
    assert(ndims==1 || ndims==2);
    H5Sget_simple_extent_dims(space_id, dims, maxdims);
    H5Sclose(space_id);
    d1 = (ndims == 1) ?       1 : dims[0];
    d2 = (ndims == 1) ? dims[0] : dims[1];
}

void read_dset_2d_d(hid_t g, const char* dname, vector<double>& v, vector<double>& w)
{
    hid_t dset_id = H5Dopen2(g, dname, H5P_DEFAULT);
    hsize_t dim1 = -1, dim2 = -1;
    get_dset2d_size(dset_id, dim1, dim2);

    v.resize(dim1); w.resize(dim1);
    hid_t memspace_id = H5Screate_simple(1, &dim1, NULL);
    hsize_t startmem[1] = {0}; hsize_t countmem[1] = {dim1};
    H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, startmem, NULL, countmem, NULL);

    hid_t filespace_id = H5Dget_space(dset_id);
    hsize_t countfile[2] = {dim1, 1}; hsize_t startfile[2] = {0, 0};
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, startfile, NULL, countfile, NULL); // Get X, mask Y
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT, &v[0]);
    startfile[1] = 1;
    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, startfile, NULL, countfile, NULL); // Get Y, mask X
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT, &w[0]);

    H5Sclose(memspace_id);
    H5Sclose(filespace_id);
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
        dset_id = H5Dcreate2(file_id, dname, H5T_STD_I32LE, dspc_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
    hid_t dataset_id = H5Dcreate2(file_id, dname, H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
    hid_t dataset_id = H5Dcreate2(file_id, dname, H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
}

void read_nodes(hid_t g, vector<double>& x, vector<double>& y)
{
    read_dset_2d_d(g, "/Nodes", x, y);
}

void read_quads(hid_t g, vector<Quad*>& v, vector<double>& x, vector<double>& y, int& t)
{
    hid_t dset_id;
    hsize_t n0, n1;
    if     (H5Lexists(g, "/Sem2D/Quad4", H5P_DEFAULT)) {dset_id=H5Dopen2(g, "/Sem2D/Quad4", H5P_DEFAULT); t=4;}
    else if(H5Lexists(g, "/Sem2D/Quad8", H5P_DEFAULT)) {dset_id=H5Dopen2(g, "/Sem2D/Quad8", H5P_DEFAULT); t=8;}
    else   {printf("ERR: only Quad4 and Quad8 are supported\n"); exit(1);}
    get_dset2d_size(dset_id, n0, n1);
    assert(t == n1);
    int* quads = new int[n0*n1];
    H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, quads);
    for(int i=0; i<n0; i++)
    {
        if (t == 4) {Quad4* q = new Quad4(quads+i*t); q->check_orient(x, y); v.push_back(q);};
        if (t == 8) {Quad8* q = new Quad8(quads+i*t); q->check_orient(x, y); v.push_back(q);};
    }
    if(quads){delete[] quads; quads=NULL;}
    H5Dclose(dset_id);
}

void read_mat(hid_t g, vector<int>& m)
{
    read_dset_1d_i(g, "/Sem2D/Mat", m);
}

void Mesh2D::read_mesh(const string& fname, int& t)
{
    hid_t file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    read_nodes(file_id, m_px, m_py);
    read_quads(file_id, m_quads, m_px, m_py, t);
    read_mat(file_id, m_mat1);
    read_mat(file_id, m_mat2);
    assert(m_mat1.size()==m_mat2.size());
    assert(m_quads.size()==m_mat1.size());
    printf("%ld Nodes, %ld Quads%i\n", m_px.size(), m_quads.size(), t);
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
        for(int k=0;k<m_quads[i]->get_nb_nodes();++k) {
            if (m_quads[i]->is_intermediate_node(k)) continue; // For topology construction, rely on principal nodes only
            eind[4*i+k] = m_quads[i]->get_node_id(k);
        }
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
    METIS_PartGraphRecursive(&ne, &ncon, xadj, adjncy, &vwgt[0], &vsize[0], &adjwgt[0], &inproc, &tpwgts[0], ubvec, options, &edgecut, &m_procs[0]);
    //    for(int k=0;k<ne;++k) printf("%d : %d\n", k, m_procs[k]);
}

void Mesh2D::write_proc_field(const string& fname)
{
    hid_t file_id = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    write_dset_1d_i(file_id, "Proc", m_procs);
    H5Fclose(file_id);
}


void Mesh2D::write_proc_file(const string& fname, const int t, int rk)
{
    hid_t fid;
    if (access(fname.c_str(), F_OK)==0) {
        // File already exists
        unlink(fname.c_str());
    }
    fid = H5Fopen(fname.c_str(), H5F_ACC_RDWR|H5F_ACC_CREAT, H5P_DEFAULT);

    MeshProcInfo info(t);
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
    write_dset_2d_i(fid, "elements", info.m_npq, info.m_quadnodes);
    write_dset_2d_i(fid, "edges", 4, info.m_edges);
    // With linear Quad, vertices==elements
    write_dset_2d_i(fid, "vertices", 4, info.m_quadvertices);
    write_dset_2d_i(fid, "faces_elem", 2, info.m_edges_elems);
    write_dset_2d_i(fid, "faces_which", 2, info.m_edges_local);
    write_dset_2d_i(fid, "faces_vertex", 2, info.m_edges_vertices);
    vector<int> vgn;
    for(map<int,int>::const_iterator it = info.m_vert_map.begin(); it != info.m_vert_map.end(); it++) {
        vgn.push_back(it->first);
    }
    write_dset_1d_i(fid, "vertices_globnum", vgn);

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

    int t = -1; // Element type : 4 or 8
    Mesh2D mesh;
    mesh.read_mesh(fmesh, t);
    if (nproc>1) {
        if (t == 8) {
            printf("Error : not yet implemented\n");
            exit(1);
        }
        mesh.partition_metis(nproc);
    }
    mesh.write_proc_field(fmesh);

    for(int k=0;k<nproc;++k) {
        snprintf(fname, 1024,"%s.%04d.h5", argv[3], k);
        mesh.write_proc_file(fname, t, k);
    }
}
// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=4 et tw=80 smartindent :*/
