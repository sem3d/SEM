/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <utility> // pair
#include <algorithm> // find
#include <tuple>
#define OMPI_SKIP_MPICXX
#include <hdf5.h>
#include <cassert>
#include <unistd.h>
#include <cstdint>
#include "metis.h"
#include "h5helper.h"
#include "read_unv.hpp"
#include <sys/stat.h>
#include <sys/types.h>
#include "../../COMMON/read_input.h"
#include "../../COMMON/sem_input.h"

using namespace std;

static void ensure_parent_dir(const string& path) {
    size_t pos = path.find_last_of('/');
    if (pos != string::npos) {
        string dir = path.substr(0, pos);
        mkdir(dir.c_str(), 0777);
    }
}

class Point {
    public:
        Point() {}
        Point(double u, double v):x(u),y(v) {}
        Point(const Point& p):x(p.x),y(p.y) {}
        double x,y;
};

typedef pair<int,int> edge_idx_t;               // Edge index : <node0,node1>                 with node0<node1
typedef tuple<int,int,int,int,int> edge_info_t; // Edge info  : <edgenum,elem0,elem1,wf0,wf1> with elemX = element X, wfX = which_face of elemX

class Quad {
    public:
        Quad() {nn = 0; n = NULL;}
        Quad(vector<uint64_t> & nodeID) {nn = nodeID.size(); n = new int[nn]; for(int i=0; i<nn; i++) n[i]=nodeID[i];}
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
            if ((v01[0]*v03[1]-v01[1]*v03[0])<0.) swap_orient();
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
        Quad4(vector<uint64_t> & nodeID):Quad(nodeID) {}
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
        Quad8(vector<uint64_t> & nodeID):Quad4(nodeID) {}
        virtual void swap_orient() {int tmp=n[4]; n[4]=n[7]; n[7]=tmp; tmp=n[5]; n[5]=n[6]; n[6]=tmp; Quad4::swap_orient();};
        virtual edge_idx_t get_edge_from_node(int i) {assert (i>=0 && i<=7); int ii = (i >= 4) ? i - 4 : i; return Quad4::get_edge_from_node(ii);};
        virtual int get_face_from_node(int i) {assert (i>=0 && i<=7); int ii = (i >= 4) ? i - 4 : i; return ii;};
        virtual bool is_intermediate_node(int i) {assert (i>=0 && i<=7); return (i >= 4) ? true : false;};
};

struct edge_comm_info_t
{
    edge_comm_info_t():edge(-1),coherency(-1) {}
    edge_comm_info_t(int e, int c):edge(e),coherency(c) {}
    edge_comm_info_t(const edge_comm_info_t& ei):edge(ei.edge),coherency(ei.coherency) {}
    int edge;
    int coherency;
};

struct Comm_proc {
    vector<int> m_vertices;
    map<int,int> m_vertices_map;
    vector<int> m_edges;
    map<edge_idx_t,edge_comm_info_t> m_edges_map;
    vector<int> m_coherency;
};

class MeshProcInfo {
    public:
        MeshProcInfo(const int npq) { m_npq = npq; }
        // Computed for storage
        vector<double> m_nodes;
        vector<int> m_quadnodes; // 4 or 8 node number for each quad
        vector<int> m_quadvertices; // always 4 vertice number for each quad
        vector<int> m_material;
        vector<int> m_quadedges; // 4 edge number for each quad

        int n_elements() const { return m_quadnodes.size()/m_npq; }
        int n_edges() const { return m_edge_map.size(); }
        int n_vertices() const { return m_vert_map.size(); }

        // Bookkeeping info
        vector<int> m_node_map; // map global node with local node (processor submesh)
        map<int, int> m_vert_map; // map node with vertex (for quad8, vertex numbering != node numbering as intermediate nodes are not considered)
        vector<int> m_quad_map;
        map<edge_idx_t, edge_info_t> m_edge_map; // maps each edge with informations associated to this edge
        vector<edge_idx_t> m_edges; // Keep track of face creation order (lost when using only map)

        // Methods
        int add_edge_element(edge_idx_t e, int fn, int qn);
        void add_vertex(const int nid) {if(m_vert_map.find(nid)==m_vert_map.end()) {int vid=n_vertices(); m_vert_map[nid]=vid;}}
        void get_quad_edge(Quad& q, int fn, int& v0, int& v1);

        map<int,Comm_proc> m_comm;
        int m_npq; // Nb nodes per quad
};

class Mesh2D {
public:
    Mesh2D() {};
    virtual ~Mesh2D() {for(vector<Quad*>::iterator it = m_quads.begin(); it != m_quads.end(); ++it) {Quad* q=*it; if(q){delete q; q=NULL;}}};
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
    vector<Quad*> m_quads;
    vector<idx_t> m_procs;
    vector<int>  m_mat1;
    vector<int>  m_mat2;
private:
    void read_sem_mesh(const string& fname);
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

int MeshProcInfo::add_edge_element(edge_idx_t e, int fn, int qn)
{
    edge_info_t ei;
    if (m_edge_map.find(e) != m_edge_map.end()) ei = m_edge_map[e];
    else {ei = edge_info_t (m_edge_map.size (), -1, -1, -1, -1); m_edges.push_back(e);}

    int elem0 = get<1>(ei); int elem1 = get<2>(ei);
    if      (elem0 == -1) {get<1>(ei) = qn; get<3>(ei) = fn;} // Modify info
    else if (elem1 == -1) {get<2>(ei) = qn; get<4>(ei) = fn;} // Modify info
    else                  {assert(0); /*should never happen - assert to make sure*/ }

    m_edge_map[e] = ei; // Store (created or modified) info in map
    return get<0>(ei);
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
            int edge_num = info.add_edge_element(quad->get_edge_from_node(k), quad->get_face_from_node(k), qn);
            info.m_quadedges.push_back(edge_num);
        }
    }
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
                auto edge_it = info.m_edge_map.find(e);
                if (edge_it != info.m_edge_map.end()) {
                    int edge_num = get<0>(edge_it->second);
                    Comm_proc& comm = info.m_comm[m_procs[qn]];
                    comm.m_edges_map[e] = edge_comm_info_t(edge_num, 1);
                }
            }
        }
    }
    map<int,Comm_proc>::iterator cit;
    for(cit=info.m_comm.begin();cit!=info.m_comm.end();++cit) {
        Comm_proc& comm=cit->second;
        for(auto eit=comm.m_edges_map.begin();eit!=comm.m_edges_map.end();++eit) {
            comm.m_edges.push_back(eit->second.edge);
            comm.m_coherency.push_back(eit->second.coherency);
        }
        map<int,int>::const_iterator vit;
        for(vit=comm.m_vertices_map.begin();vit!=comm.m_vertices_map.end();++vit) {
            comm.m_vertices.push_back(vit->second);
        }
    }
}


void read_nodes(hid_t g, vector<double>& x, vector<double>& y)
{
    h5h_read_dset_Nx2(g, "/Nodes", x, y);
}

void read_quads(hid_t g, vector<Quad*>& v, vector<double>& x, vector<double>& y)
{
    hid_t dset_id;
    hsize_t n0, n1;
    int t = 0;
    if     (H5Lexists(g, "/Sem2D/Quad4", H5P_DEFAULT)) {dset_id=H5Dopen2(g, "/Sem2D/Quad4", H5P_DEFAULT); t=4;}
    else if(H5Lexists(g, "/Sem2D/Quad8", H5P_DEFAULT)) {dset_id=H5Dopen2(g, "/Sem2D/Quad8", H5P_DEFAULT); t=8;}
    else   {printf("ERR: only Quad4 and Quad8 are supported\n"); exit(1);}
    h5h_get_dset2d_size(dset_id, n0, n1);
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
    h5h_read_dset(g, "/Sem2D/Mat", m);
}

void Mesh2D::read_mesh(const string& fname)
{
  if (fname.find(".unv") != string::npos)
  {
    vector<int> filterelems; filterelems.push_back(44); filterelems.push_back(45);
    lsnodes nodes; lselems elems;
    int rc = read_unv_mesh(fname, nodes, elems, &filterelems);
    assert(rc == 0);

    // Look for nodes in the XY plan (UNV, which is a 3D generic format, do not provide any hint to "choose" between XY or XZ or YZ)

    for (unsigned int i = 0; i < nodes.size(); i++) {m_px.push_back(get<0>(nodes[i])); m_py.push_back(get<1>(nodes[i]));}

    // Look for quad4 / quad8

    for (unsigned int i = 0; i < elems.size(); i++)
    {
      int type = get<0>(elems[i]);
      vector<uint64_t> nodeID = get<2>(elems[i]);
      group gp = get<3>(elems[i]);
      if(type == 44) {Quad4* q = new Quad4(nodeID); q->check_orient(m_px, m_py); m_quads.push_back(q); m_mat1.push_back(get<1>(gp));}; // Quad4
      if(type == 45) {Quad8* q = new Quad8(nodeID); q->check_orient(m_px, m_py); m_quads.push_back(q); m_mat1.push_back(get<1>(gp));}; // Quad8
    }
    m_mat2 = m_mat1;
  }
  else read_sem_mesh(fname);

  m_nprocs = 1;
  m_procs.resize(m_quads.size(), 0);
  m_mat_max = 0;
  for(int k=0;k<m_mat1.size();++k) if (m_mat1[k]>m_mat_max) m_mat_max = m_mat1[k];

  assert(m_px.size() > 0 && m_px.size() == m_py.size());               // Check nodes    consistency
  assert(m_quads.size() > 0);                                          // Check element  consistency
  assert(m_mat1.size() > 0 && m_mat1.size() == m_mat2.size());         // Check material consistency
  assert(std::find(m_mat1.begin(), m_mat1.end(), -1) == m_mat1.end()); // Check all elements have a material
}

void Mesh2D::read_sem_mesh(const string& fname)
{
    hid_t file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    read_nodes(file_id, m_px, m_py);
    read_quads(file_id, m_quads, m_px, m_py);
    read_mat(file_id, m_mat1);
    read_mat(file_id, m_mat2);
    assert(m_mat1.size()==m_mat2.size());
    assert(m_quads.size()==m_mat1.size());
    assert(m_quads.size()>0);
    printf("%ld Nodes, %ld Quads%i\n", m_px.size(), m_quads.size(), m_quads[0]->get_nb_nodes());
    H5Fclose(file_id);
}

void Mesh2D::check_cell_orient()
{
}

void Mesh2D::partition_metis(int nproc)
{
    if (m_quads.size() > 0 && m_quads[0]->get_nb_nodes() == 8) {
        printf("Error : not yet implemented\n");
        exit(1);
    }

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
    hid_t file_id = H5Fopen(fname.c_str(), H5F_ACC_RDWR|H5F_ACC_CREAT, H5P_DEFAULT);
    h5h_write_dset(file_id, "Proc", m_procs);
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

    assert(m_quads.size() > 0);
    MeshProcInfo info(m_quads[0]->get_nb_nodes());
    gather_proc_info(info, rk);
    h5h_write_attr_int(fid, "ndim", 2);
    //m_nprocs=1;
    h5h_write_attr_int(fid, "n_processors", m_nprocs);
    h5h_write_attr_int(fid, "n_materials", m_mat_max+1);
    h5h_write_attr_int(fid, "n_elements", info.n_elements());
    h5h_write_attr_int(fid, "n_edges", info.n_edges());
    h5h_write_attr_int(fid, "n_vertices", info.n_vertices());
    h5h_write_dset_2d(fid, "nodes", 2, info.m_nodes);
    h5h_write_dset_2d(fid, "material", 3, info.m_material);
    h5h_write_dset_2d(fid, "elements", info.m_npq, info.m_quadnodes);
    h5h_write_dset_2d(fid, "edges", 4, info.m_quadedges);
    // With linear Quad, vertices==elements
    h5h_write_dset_2d(fid, "vertices", 4, info.m_quadvertices);
    map<int, int> elem_map;
    int local_elem_counter = 0;
    for(int qn=0; qn<m_quads.size(); ++qn) {
        if (rk == -1 || rk == m_procs[qn]) {
            elem_map[qn] = local_elem_counter++;
        }
    }

    vector<int> edges_elems, edges_wf, edges_vertices;
    for(auto it = info.m_edges.begin(); it != info.m_edges.end(); it++) {
        edge_idx_t e = *it;
        edges_vertices.push_back(info.m_vert_map[info.m_node_map[get<0>(e)]]);
        edges_vertices.push_back(info.m_vert_map[info.m_node_map[get<1>(e)]]);
        edge_info_t ei = info.m_edge_map[e];
        int el0 = get<1>(ei);
        int el1 = get<2>(ei);
        if (el0 != -1 && elem_map.find(el0) != elem_map.end()) {
            el0 = elem_map[el0];
        } else {
            el0 = -1;
        }
        if (el1 != -1 && elem_map.find(el1) != elem_map.end()) {
            el1 = elem_map[el1];
        } else {
            el1 = -1;
        }
        edges_elems.push_back(el0);
        edges_elems.push_back(el1);
        edges_wf.push_back(get<3>(ei));
        edges_wf.push_back(get<4>(ei));
    }
    h5h_write_dset_2d(fid, "faces_elem", 2, edges_elems);
    h5h_write_dset_2d(fid, "faces_which", 2, edges_wf);
    h5h_write_dset_2d(fid, "faces_vertex", 2, edges_vertices);
    vector<int> vgn;
    for(map<int,int>::const_iterator it = info.m_vert_map.begin(); it != info.m_vert_map.end(); it++) {
        vgn.push_back(it->first);
    }
    h5h_write_dset(fid, "vertices_globnum", vgn);

    int n_comm = info.m_comm.size();
    //n_comm = 0;
    h5h_write_attr_int(fid, "n_communications", n_comm);
    map<int,Comm_proc>::const_iterator it;
    int comm_count=0;
    for(it=info.m_comm.begin();it!=info.m_comm.end();++it) {
        char grp_name[60];
        snprintf(grp_name, 60, "Comm%05d", comm_count);
        const Comm_proc& comm = it->second;
        hid_t grp = H5Gcreate(fid, grp_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        h5h_write_attr_int(grp, "processor", it->first);
        h5h_write_dset_empty(grp, "vertices", comm.m_vertices);
        h5h_write_dset_empty(grp, "edges", comm.m_edges);
        h5h_write_dset_empty(grp, "coherency", comm.m_coherency);
        H5Gclose(grp);
        comm_count++;
    }

    H5Fclose(fid);

    if (info.n_elements() <= 0) {
        return;
    }

    // Automatically generate a companion .xmf file for visual inspection of the mesh in ParaView
    string xmf_name = fname;
    size_t h5_pos = xmf_name.rfind(".h5");
    if (h5_pos != string::npos) {
        xmf_name.replace(h5_pos, 3, ".xmf");
    } else {
        xmf_name += ".xmf";
    }

    string h5_filename = fname;
    size_t slash_pos = h5_filename.find_last_of('/');
    if (slash_pos != string::npos) {
        h5_filename = h5_filename.substr(slash_pos + 1);
    }

    FILE* fx = fopen(xmf_name.c_str(), "w");
    if (fx) {
        fprintf(fx, "<?xml version=\"1.0\" ?>\n");
        fprintf(fx, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n");
        fprintf(fx, "<Xdmf Version=\"2.0\">\n");
        fprintf(fx, "  <Domain>\n");
        fprintf(fx, "    <Grid Name=\"mesh.%04d\" GridType=\"Uniform\">\n", rk);
        fprintf(fx, "      <Topology Type=\"Quadrilateral\" NumberOfElements=\"%d\">\n", info.n_elements());
        fprintf(fx, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%d 4\">\n", info.n_elements());
        fprintf(fx, "          %s:/elements\n", h5_filename.c_str());
        fprintf(fx, "        </DataItem>\n");
        fprintf(fx, "      </Topology>\n");
        fprintf(fx, "      <Geometry Type=\"XY\">\n");
        fprintf(fx, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%d 2\">\n", info.n_vertices());
        fprintf(fx, "          %s:/nodes\n", h5_filename.c_str());
        fprintf(fx, "        </DataItem>\n");
        fprintf(fx, "      </Geometry>\n");
        fprintf(fx, "      <Attribute Name=\"MaterialID\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(fx, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%d 3\" ItemType=\"HyperSlab\">\n", info.n_elements());
        fprintf(fx, "          <DataItem Format=\"XML\" Dimensions=\"3 3\">\n");
        fprintf(fx, "            0 0 0\n");
        fprintf(fx, "            1 1 1\n");
        fprintf(fx, "            %d 1 3\n", info.n_elements());
        fprintf(fx, "          </DataItem>\n");
        fprintf(fx, "          <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%d 3\">\n", info.n_elements());
        fprintf(fx, "            %s:/material\n", h5_filename.c_str());
        fprintf(fx, "          </DataItem>\n");
        fprintf(fx, "        </DataItem>\n");
        fprintf(fx, "      </Attribute>\n");
        fprintf(fx, "    </Grid>\n");
        fprintf(fx, "  </Domain>\n");
        fprintf(fx, "</Xdmf>\n");
        fclose(fx);
    }
}

// Read a line from f, skipping comment lines that start with '#'.
// Mirrors the 3D mesher's getData_line so 2D mat.dat files share the syntax.
static void getData_line(char** buffer, size_t* linesize, FILE* f)
{
    ssize_t nc = getline(buffer, linesize, f);
    for(int k=0;k<100;++k) {
        if((*buffer)[0] != '#') break;
        nc = getline(buffer, linesize, f);
    }
    if (nc<=0 && *buffer) (*buffer)[0] = 0;
}

// One SEM2D material as carried by mater.in / material.input.
// The 2D solver (Domain.F90 read_material_file) needs NGLL and Dt per material, unlike the
// 3D mesher which reads those from input.spec -- so the 2D material file format differs.
struct Material2D {
    char   type;            // S solid, F fluid, P pml
    double vp, vs, rho;
    int    ngll;
    double dt;
    double qp, qs;
    // PML descriptor (type=='P' only), as read by Domain.F90 read_material_file:
    //   Filtering npow Apow Px Left Pz Down omegac kc
    bool   is_pml;
    bool   px, left, pz, down;  // PML in x / on x-min side ; PML in z / on z-min side
    int    npow;
    double apow, omegac, kc;
    Material2D() : type('S'), vp(0), vs(0), rho(0), ngll(5), dt(0), qp(0), qs(0),
                   is_pml(false), px(false), left(false), pz(false), down(false),
                   npow(2), apow(10.), omegac(0.), kc(0.) {}
};

// mater.in (2D): a header-commented file, then
//   n_mat
//   <type> <Vp> <Vs> <Rho> <NGLL> <Dt> <Qp> <Qs>      (n_mat lines)
static void read_materials_2d(const char* fname, vector<Material2D>& mats)
{
    FILE* f = fopen(fname, "r");
    if (!f) {printf("ERR: cannot open %s\n", fname); exit(1);}
    char* buffer=NULL; size_t n=0;
    int nmat=0;
    getData_line(&buffer, &n, f);
    sscanf(buffer, "%d", &nmat);
    if (nmat<=0 || nmat>1000) {printf("ERR: bad material count %d in %s\n", nmat, fname); exit(1);}
    for(int k=0;k<nmat;++k) {
        Material2D m;
        m.type='S'; m.vp=m.vs=m.rho=0.; m.ngll=5; m.dt=0.; m.qp=0.; m.qs=0.;
        getData_line(&buffer, &n, f);
        int c = sscanf(buffer, " %c %lf %lf %lf %d %lf %lf %lf",
                       &m.type, &m.vp, &m.vs, &m.rho, &m.ngll, &m.dt, &m.qp, &m.qs);
        if (c<5) {printf("ERR: material line %d in %s has too few fields (%d)\n", k, fname, c); exit(1);}
        mats.push_back(m);
    }
    if(buffer) free(buffer);
    fclose(f);
}

// material.input (2D), read by Domain.F90 read_material_file:
//   n_mat
//   <type> <Vp> <Vs> <Rho> <NGLLx> <mid> <NGLLz> <Dt> <Qp> <Qs>      (n_mat lines)
//   [if any 'P' material:]
//   <two comment/header lines>
//   <Filtering npow Apow Px Left Pz Down omegac kc>                  (one per 'P', in order)
static void write_materials_2d(const char* fname, const vector<Material2D>& mats)
{
    FILE* f = fopen(fname, "w");
    if (!f) {printf("ERR: cannot write %s\n", fname); exit(1);}
    fprintf(f, "%ld\n", mats.size());
    int npml = 0;
    for(size_t k=0;k<mats.size();++k) {
        const Material2D& m = mats[k];
        fprintf(f, "%c %g %g %g %d %d %d %g %g %g\n",
                m.type, m.vp, m.vs, m.rho,
                m.ngll, m.ngll, m.ngll, m.dt, m.qp, m.qs);
        if (m.is_pml) npml++;
    }
    if (npml>0) {
        fprintf(f, "# Specifications for PMLs\n");
        fprintf(f, "# Filtering, npow, Apow, Px, Left, Pz, Down, omegac, kc\n");
        for(size_t k=0;k<mats.size();++k) {
            const Material2D& m = mats[k];
            if (!m.is_pml) continue;
            fprintf(f, "F %d %g %c %c %c %c %g %g\n",
                    m.npow, m.apow,
                    m.px?'T':'F', m.left?'T':'F', m.pz?'T':'F', m.down?'T':'F',
                    m.omegac, m.kc);
        }
    }
    fclose(f);
}

// Axis-aligned 2D rectangular grid built from mat.dat (the 2D analog of the 3D RectMesh).
// The vertical (layered) direction is z, stored in Mesh2D::m_py.
struct RectMesh2D {
    double xmin, xmax, xstep, zmax;
    int    nlayers;
    vector<double> thickness;
    vector<int>    nsteps;
    int    npml;          // number of PML element layers (0 = no PML)
    bool   pml_W, pml_E, pml_U, pml_D; // PML present on each side
    int    ngll_pml;      // NGLL for PML elements (<=0 -> use base material NGLL)
    int    npow;          // PML attenuation exponent
    double apow, omegac, kc;
    int    elem_shape;    // 4 = Quad4
    int    nelemx, nelemz;
    // PML material cache: key (layer,W,E,U,D) -> material index (appended to mats)
    map<int,int> pml_cache;

    void read_params(const char* fname);
    void apply_pml_borders();
    int  get_mat(vector<Material2D>& mats, int layer, bool W, bool E, bool U, bool D);
    void build(Mesh2D& mesh, vector<Material2D>& mats);
    int  pointidx(int i, int zlev) const { return i + zlev*(nelemx+1); }
};

// mat.dat (2D) = the 3D mat.dat without the y-block:
//   xmin / xmax / xstep / zmax / nlayers / (thickness nsteps) x nlayers /
//   has_pml (npml) / pml_top pml_bottom / ngllPML [npow Apow omegac kc] / mesh_type
// Lateral PML (W,E) is always on when npml>0 (as in the 3D mesher default). pml_top/pml_bottom
// toggle the U/D sides. When has_pml=0 the next two lines are dummies.
void RectMesh2D::read_params(const char* fname)
{
    FILE* f = fopen(fname, "r");
    if (!f) {printf("ERR: cannot open %s\n", fname); exit(1);}
    char* buffer=NULL; size_t n=0;

    getData_line(&buffer,&n,f); sscanf(buffer, "%lf", &xmin);
    getData_line(&buffer,&n,f); sscanf(buffer, "%lf", &xmax);
    getData_line(&buffer,&n,f); sscanf(buffer, "%lf", &xstep);
    getData_line(&buffer,&n,f); sscanf(buffer, "%lf", &zmax);
    getData_line(&buffer,&n,f); sscanf(buffer, "%d",  &nlayers);
    if (nlayers<1 || nlayers>200) {printf("ERR: bad nlayers=%d in %s\n", nlayers, fname); exit(1);}
    thickness.resize(nlayers); nsteps.resize(nlayers);
    for(int k=0;k<nlayers;++k) {
        getData_line(&buffer,&n,f);
        sscanf(buffer, "%lf %d", &thickness[k], &nsteps[k]);
        if (thickness[k]<=0. || nsteps[k]<1) {printf("ERR: bad layer %d (thick=%g nsteps=%d)\n", k, thickness[k], nsteps[k]); exit(1);}
    }

    npml = 0;
    getData_line(&buffer,&n,f); sscanf(buffer, "%d", &npml);
    if (npml<0) {printf("ERR: has_pml=%d must be >=0\n", npml); exit(1);}

    int pml_top=0, pml_bottom=1;
    getData_line(&buffer,&n,f); sscanf(buffer, "%d %d", &pml_top, &pml_bottom);

    ngll_pml=0; npow=2; apow=10.; omegac=0.; kc=0.;
    getData_line(&buffer,&n,f); sscanf(buffer, "%d %d %lf %lf %lf", &ngll_pml, &npow, &apow, &omegac, &kc);

    getData_line(&buffer,&n,f); elem_shape=4; sscanf(buffer, "%d", &elem_shape);
    if (elem_shape!=4) {printf("ERR: only mesh_type 4 (Quad4) is supported on the fly\n"); exit(1);}

    pml_W = pml_E = (npml>0);
    pml_U = (npml>0 && pml_top);
    pml_D = (npml>0 && pml_bottom);

    if(buffer) free(buffer);
    fclose(f);
}

// Extend the domain outward by npml element layers on each active PML side (mirrors the 3D
// RectMesh::apply_pml_borders). The top/bottom extension grows the first/last layer.
void RectMesh2D::apply_pml_borders()
{
    if (npml<=0) return;
    if (pml_E) xmax += npml*xstep;
    if (pml_W) xmin -= npml*xstep;
    if (pml_U) {
        double zstep = thickness[0]/nsteps[0];
        zmax        += npml*zstep;
        thickness[0]+= npml*zstep;
        nsteps[0]   += npml;
    }
    if (pml_D) {
        int ll = nlayers-1;
        double zstep = thickness[ll]/nsteps[ll];
        thickness[ll]+= npml*zstep;
        nsteps[ll]   += npml;
    }
}

// Return the material index for an element in layer `layer` touching the given PML sides.
// Non-PML elements keep their layer index; PML elements get a derived 'P' material (created
// once per (layer, side-combination) and cached).
int RectMesh2D::get_mat(vector<Material2D>& mats, int layer, bool W, bool E, bool U, bool D)
{
    bool px = W || E;
    bool pz = U || D;
    if (!px && !pz) return layer;

    int key = layer*16 + (W?1:0) + (E?2:0) + (U?4:0) + (D?8:0);
    map<int,int>::iterator it = pml_cache.find(key);
    if (it != pml_cache.end()) return it->second;

    Material2D m = mats[layer];   // copy base properties (vp/vs/rho/dt)
    m.type   = 'P';
    m.ngll   = (ngll_pml>0) ? ngll_pml : mats[layer].ngll;
    m.is_pml = true;
    m.px = px;  m.left = W;        // left side = x-min
    m.pz = pz;  m.down = D;        // down side = z-min
    m.npow = npow;  m.apow = apow;  m.omegac = omegac;  m.kc = kc;

    int idx = mats.size();
    mats.push_back(m);
    pml_cache[key] = idx;
    return idx;
}

void RectMesh2D::build(Mesh2D& mesh, vector<Material2D>& mats)
{
    apply_pml_borders();
    if (!(xmin<xmax)) {printf("ERR: need xmin<xmax (%g,%g)\n", xmin, xmax); exit(1);}

    nelemx = int((xmax-xmin)/xstep);
    nelemz = 0;
    for(int k=0;k<nlayers;++k) nelemz += nsteps[k];
    if (nelemx<1 || nelemz<1) {printf("ERR: empty grid %dx%d\n", nelemx, nelemz); exit(1);}
    printf("Creating grid mesh %d x %d with linear (Quad4) elements (npml=%d)\n", nelemx, nelemz, npml);

    // Nodes: z-levels from top (zmax) downward, shared interface rows not duplicated.
    double layerzmax = zmax;
    int k0 = 0;
    for(int nl=0;nl<nlayers;++nl) {
        double zstep = thickness[nl]/nsteps[nl];
        for(int k=k0;k<=nsteps[nl];++k) {
            double z = layerzmax - k*zstep;
            for(int i=0;i<=nelemx;++i) {
                mesh.m_px.push_back(xmin + i*xstep);
                mesh.m_py.push_back(z);
            }
        }
        k0 = 1;
        layerzmax -= thickness[nl];
    }

    // Elements: zlev increases downward. PML sides flagged on the outer npml rings.
    int k = 0; // global z-cell index
    for(int nl=0;nl<nlayers;++nl) {
        for(int kl=0;kl<nsteps[nl];++kl) {
            for(int i=0;i<nelemx;++i) {
                int q[4];
                q[0] = pointidx(i,   k+1); // bottom-left
                q[1] = pointidx(i+1, k+1); // bottom-right
                q[2] = pointidx(i+1, k  ); // top-right
                q[3] = pointidx(i,   k  ); // top-left
                Quad4* qd = new Quad4(q);
                qd->check_orient(mesh.m_px, mesh.m_py);
                mesh.m_quads.push_back(qd);

                bool W = pml_W && (i < npml);
                bool E = pml_E && (i > nelemx-npml-1);
                bool U = pml_U && (k < npml) && (nl==0);
                bool D = pml_D && (k > nelemz-npml-1) && (nl==nlayers-1);
                mesh.m_mat1.push_back(get_mat(mats, nl, W,E,U,D));
            }
            k++;
        }
    }
    mesh.m_mat2 = mesh.m_mat1;
    mesh.m_nprocs = 1;
    mesh.m_procs.assign(mesh.m_quads.size(), 0);
    mesh.m_mat_max = 0;
    for(size_t e=0;e<mesh.m_mat1.size();++e)
        if (mesh.m_mat1[e]>mesh.m_mat_max) mesh.m_mat_max = mesh.m_mat1[e];
}

// On-the-fly axis-aligned quad grid built from mat.dat + mater.in; writes material.input.
void handle_on_the_fly(Mesh2D& mesh)
{
    vector<Material2D> mats;
    read_materials_2d("mater.in", mats);

    RectMesh2D desc;
    desc.read_params("mat.dat");
    if (desc.nlayers > (int)mats.size()) {
        printf("ERR: mat.dat has %d layers but mater.in only defines %ld materials\n",
               desc.nlayers, mats.size());
        exit(1);
    }
    desc.build(mesh, mats);   // appends derived PML materials to mats

    write_materials_2d("material.input", mats);
    printf("Wrote material.input (%ld materials)\n", mats.size());
}

// ===================================================================
// PML by extrusion of boundary edges, driven by an optional pml.input.
// The 2D analog of the 3D mesher's mesh_pml_extrude. Sides: x-, x+, z-, z+
// (z is the vertical axis, stored in Mesh2D::m_py). Processed x then z so a
// W+D corner PML is created when the z- pass meets the z-min edge of an
// x-column. y-/y+ are 3D-only and rejected.
enum Pml2Side { P2_XM=0, P2_XP, P2_ZM, P2_ZP, P2_NSIDES };

struct PmlSpec2D {
    int    n[P2_NSIDES];
    double step[P2_NSIDES];
    int    npow;
    double apow, omegac, kc;
    PmlSpec2D():npow(2),apow(10.),omegac(0.),kc(0.) {
        for(int k=0;k<P2_NSIDES;++k){ n[k]=0; step[k]=0.; }
    }
    bool any() const { for(int k=0;k<P2_NSIDES;++k) if(n[k]>0) return true; return false; }
};

static bool read_pml_input_2d(const char* fname, PmlSpec2D& spec)
{
    FILE* f = fopen(fname, "r");
    if (!f) return false;
    char* buffer=NULL; size_t n=0;
    while (true) {
        getData_line(&buffer, &n, f);
        if (!buffer || buffer[0]==0) break;
        char tok[64]={0};
        if (sscanf(buffer, "%63s", tok)!=1) continue;
        if (!strcmp(tok,"pmlparams")) {
            sscanf(buffer, "%*s %d %lf %lf %lf", &spec.npow, &spec.apow, &spec.omegac, &spec.kc);
            continue;
        }
        int nn=0; double step=0.;
        int c = sscanf(buffer, "%*s %d %lf", &nn, &step);
        int s=-1;
        if (!strcmp(tok,"x-")) s=P2_XM;
        else if (!strcmp(tok,"x+")) s=P2_XP;
        else if (!strcmp(tok,"z-")) s=P2_ZM;
        else if (!strcmp(tok,"z+")) s=P2_ZP;
        else if (!strcmp(tok,"y-") || !strcmp(tok,"y+")) {
            printf("ERR pml.input: side '%s' is 3D-only; 2D has no y axis\n", tok); exit(1);
        } else { printf("ERR pml.input: unknown side '%s'\n", tok); exit(1); }
        if (c<1 || nn<0) { printf("ERR pml.input: bad count for '%s'\n", tok); exit(1); }
        spec.n[s]=nn; spec.step[s]=(c>=2)?step:0.;
    }
    if (buffer) free(buffer);
    printf("Read pml.input: x-=%d x+=%d z-=%d z+=%d\n",
           spec.n[P2_XM], spec.n[P2_XP], spec.n[P2_ZM], spec.n[P2_ZP]);
    return true;
}

struct MatInfo2D { int base; bool W,E,D,U; MatInfo2D():base(-1),W(false),E(false),D(false),U(false){} };

struct PmlExtruder2D {
    Mesh2D& mesh;
    vector<Material2D>& mats;
    const PmlSpec2D& spec;
    vector<MatInfo2D> minfo;
    map<int,int> matcache;              // (base<<4|flags) -> material index
    map<pair<int,int>,int> newnode;     // (orig node, layer) -> new node id

    PmlExtruder2D(Mesh2D& m, vector<Material2D>& mt, const PmlSpec2D& s):mesh(m),mats(mt),spec(s) {
        minfo.resize(mats.size());
        for(size_t k=0;k<mats.size();++k) {
            MatInfo2D& mi=minfo[k];
            if (mats[k].is_pml) { mi.base=(int)k; mi.W=mats[k].left&&mats[k].px; mi.E=mats[k].px&&!mats[k].left;
                                  mi.D=mats[k].down&&mats[k].pz; mi.U=mats[k].pz&&!mats[k].down; }
            else mi.base=(int)k;
        }
    }

    double coord(int node, int axis) const { return axis==0 ? mesh.m_px[node] : mesh.m_py[node]; }

    int newpt(int orig, int layer, int axis, double delta) {
        pair<int,int> key(orig,layer);
        map<pair<int,int>,int>::iterator it=newnode.find(key);
        if (it!=newnode.end()) return it->second;
        double x=mesh.m_px[orig], y=mesh.m_py[orig];
        if (axis==0) x += delta*layer; else y += delta*layer;
        int id = mesh.m_px.size();
        mesh.m_px.push_back(x); mesh.m_py.push_back(y);
        newnode[key]=id;
        return id;
    }

    int get_or_make_pml(int src_mat, int side) {
        MatInfo2D si = minfo[src_mat];
        int base=si.base;
        bool W=si.W,E=si.E,D=si.D,U=si.U;
        bool* flag[P2_NSIDES]={&W,&E,&D,&U};
        if (*flag[side]) return src_mat;
        *flag[side]=true;
        bool px=W||E, left=W, pz=D||U, down=D;
        int key = (base<<4) | (px?1:0)|(left?2:0)|(pz?4:0)|(down?8:0);
        map<int,int>::iterator it=matcache.find(key);
        if (it!=matcache.end()) return it->second;
        Material2D m = mats[base];
        m.type='P'; m.is_pml=true;
        m.px=px; m.left=left; m.pz=pz; m.down=down;
        m.npow=spec.npow; m.apow=spec.apow; m.omegac=spec.omegac; m.kc=spec.kc;
        int idx=mats.size();
        mats.push_back(m);
        MatInfo2D ni; ni.base=base; ni.W=W;ni.E=E;ni.D=D;ni.U=U;
        minfo.push_back(ni);
        matcache[key]=idx;
        printf("  PML material %d (base %d) flags W%d E%d D%d U%d\n", idx, base, W,E,D,U);
        return idx;
    }

    void emit_quad(int i0, int i1, int o0, int o1, int mat) {
        int q[4] = { i0, i1, o1, o0 };      // around the quad; check_orient fixes winding
        Quad4* qd = new Quad4(q);
        qd->check_orient(mesh.m_px, mesh.m_py);
        mesh.m_quads.push_back(qd);
        mesh.m_mat1.push_back(mat);
    }

    void extrude_side(int side) {
        int nlay = spec.n[side];
        if (nlay<=0) return;
        int axis = (side==P2_XM||side==P2_XP) ? 0 : 1;
        double sgn = (side==P2_XM||side==P2_ZM) ? -1. : 1.;
        newnode.clear();
        // bbox along axis
        double lo=coord(0,axis), hi=lo;
        size_t nn = mesh.m_px.size();
        for(size_t k=0;k<nn;++k){ double c=coord((int)k,axis); if(c<lo)lo=c; if(c>hi)hi=c; }
        double ext=hi-lo; if(ext<=0.){printf("ERR: degenerate 2D mesh\n"); exit(1);}
        double plane=(sgn<0.)?lo:hi;
        double tol=1e-6*ext;
        // Quad edges as ordered node-index pairs (local): bottom,right,top,left
        static const int EDGE[4][2]={{0,1},{1,2},{3,2},{0,3}};
        struct BEdge { int a,b,mat; };
        vector<BEdge> edges;
        double step=spec.step[side], tstep=0.;
        size_t nq0 = mesh.m_quads.size();
        for(size_t e=0;e<nq0;++e) {
            Quad* q=mesh.m_quads[e];
            if (q->get_nb_nodes()!=4) { printf("ERR: 2D PML extrusion supports only Quad4\n"); exit(1); }
            // element extent along axis
            double emin=coord(q->get_node_id(0),axis), emax=emin;
            for(int k=1;k<4;++k){ double c=coord(q->get_node_id(k),axis); if(c<emin)emin=c; if(c>emax)emax=c; }
            for(int ed=0;ed<4;++ed) {
                int a=q->get_node_id(EDGE[ed][0]), b=q->get_node_id(EDGE[ed][1]);
                if (fabs(coord(a,axis)-plane)<=tol && fabs(coord(b,axis)-plane)<=tol) {
                    double th=emax-emin;
                    if (tstep<=0.) tstep=th;
                    else if (step<=0. && fabs(th-tstep)>1e-3*tstep) {
                        printf("ERR: boundary elements on side %d non-uniform (%g vs %g); set explicit step\n",
                               side, th, tstep); exit(1);
                    }
                    BEdge be={a,b,mesh.m_mat1[e]}; edges.push_back(be);
                }
            }
        }
        if (edges.empty()) { printf("WARNING: no boundary edges on side %d, skipped\n", side); return; }
        if (step<=0.) step=tstep;
        double delta=sgn*step;
        printf("Side %d: %zu boundary edges, %d layers, step=%g\n", side, edges.size(), nlay, step);
        for(size_t i=0;i<edges.size();++i) {
            int mat = get_or_make_pml(edges[i].mat, side);
            int a=edges[i].a, b=edges[i].b;
            for(int l=1;l<=nlay;++l) {
                int i0=(l==1)?a:newpt(a,l-1,axis,delta);
                int i1=(l==1)?b:newpt(b,l-1,axis,delta);
                int o0=newpt(a,l,axis,delta);
                int o1=newpt(b,l,axis,delta);
                emit_quad(i0,i1,o0,o1,mat);
            }
        }
    }

    void run() { for(int s=0;s<P2_NSIDES;++s) extrude_side(s); }
};

// Read mater.in as the base materials, add PML by extrusion per pml.input,
// rewrite material.input. Returns true if extrusion happened.
static bool extrude_pml_from_input(Mesh2D& mesh)
{
    PmlSpec2D spec;
    if (!read_pml_input_2d("pml.input", spec) || !spec.any()) return false;
    vector<Material2D> mats;
    read_materials_2d("mater.in", mats);
    if (mesh.m_mat_max >= (int)mats.size()) {
        printf("ERR: mesh references material %d but mater.in defines only %ld\n",
               mesh.m_mat_max, mats.size());
        exit(1);
    }
    printf("Extruding 2D PML layers onto imported mesh...\n");
    PmlExtruder2D ex(mesh, mats, spec);
    ex.run();
    // Refresh derived mesh bookkeeping after adding quads.
    mesh.m_mat2 = mesh.m_mat1;
    mesh.m_procs.assign(mesh.m_quads.size(), 0);
    mesh.m_mat_max = 0;
    for(size_t e=0;e<mesh.m_mat1.size();++e)
        if (mesh.m_mat1[e]>mesh.m_mat_max) mesh.m_mat_max=mesh.m_mat1[e];
    write_materials_2d("material.input", mats);
    printf("Wrote material.input (%ld materials)\n", mats.size());
    return true;
}

// Mesh2D::read_mesh ingests a single file (its own node/element numbering); it is not a
// merging reader, so we read exactly one .unv / HDF5 file (as the original mesh2dc did).
void handle_ideas_file(Mesh2D& mesh)
{
    char fname[2048];
    printf("\n.unv file name ?\n");
    if (scanf("%2000s", fname)!=1) {printf("ERR: cannot read file name\n"); exit(1);}
    mesh.read_mesh(fname);
}

void handle_hdf5_file(Mesh2D& mesh)
{
    char fname[2048];
    printf("\nHDF5 mesh file name ?\n");
    if (scanf("%2000s", fname)!=1) {printf("ERR: cannot read file name\n"); exit(1);}
    mesh.read_mesh(fname);
}

// Emulates the 3D mesher interface: stdin-driven menu, fixed output base name.
int main(int argc, char** argv)
{
    char fname[1024];
    int NPROCS;
    int choice;
    char* buffer=NULL;
    size_t linesize=0;
    string basename = "mesh4spec"; // matches mesh_file in input.spec
    if (access("input.spec", F_OK) == 0) {
        sem_config_t config;
        int err;
        read_sem_config(&config, 0, 2, "input.spec", &err);
        if (err > 0 && config.mesh_file) {
            basename = config.mesh_file;
        }
    }

    printf("-------------------------------------------------\n");
    printf("-----                                       -----\n");
    printf("----- Construction of input files for SEM2D -----\n");
    printf("-----                                       -----\n");
    printf("-------------------------------------------------\n");

    printf("\n   --> How many procs for the run ?\n");
    getData_line(&buffer, &linesize, stdin);
    sscanf(buffer, "%d", &NPROCS);
    printf("             %d processor(s)\n\n", NPROCS);

    printf("  --> Which Initial Mesh?\n");
    printf("      1- On the fly\n");
    printf("      3- Ideas (.unv) files\n");
    printf("      4- HDF5 Quad files\n");
    getData_line(&buffer, &linesize, stdin);
    sscanf(buffer, "%d", &choice);
    printf("            Your choice is %d \n\n", choice);

    Mesh2D mesh;
    switch(choice) {
    case 1:
        handle_on_the_fly(mesh);
        break;
    case 3:
        handle_ideas_file(mesh);
        break;
    case 4:
        handle_hdf5_file(mesh);
        break;
    default:
        printf("ERR: unknown mesh source %d\n", choice);
        exit(1);
    };

    // Optional PML extrusion for imported meshes (cases 3/4). Case 1 already
    // builds its PML from mat.dat, so pml.input is ignored there.
    if (choice==3 || choice==4) {
        extrude_pml_from_input(mesh);
    } else if (choice==1 && access("pml.input", F_OK)==0) {
        printf("WARNING: pml.input ignored for 'on the fly' meshes (PML comes from mat.dat)\n");
    }

    if (NPROCS>1) mesh.partition_metis(NPROCS);

    ensure_parent_dir(basename);
    mesh.write_proc_field(basename+".h5");
    for(int k=0;k<NPROCS;++k) {
        snprintf(fname, 1024,"%s.%04d.h5", basename.c_str(), k);
        mesh.write_proc_file(fname, k);
    }
    if(buffer){free(buffer); buffer=NULL;}
    return 0;
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
