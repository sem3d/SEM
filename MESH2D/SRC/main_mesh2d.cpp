/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include <cstdio>
#include <cstdlib>
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

using namespace std;

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
                int edge_num = get<0>(info.m_edge_map[e]);
                if (edge_num>0) {
                    Comm_proc& comm = info.m_comm[m_procs[qn]];
                    comm.m_edges_map[e]=edge_comm_info_t(edge_num,1);
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
    vector<int> edges_elems, edges_wf, edges_vertices;
    for(auto it = info.m_edges.begin(); it != info.m_edges.end(); it++) {
        edge_idx_t e = *it;
        edges_vertices.push_back(info.m_vert_map[info.m_node_map[get<0>(e)]]);
        edges_vertices.push_back(info.m_vert_map[info.m_node_map[get<1>(e)]]);
        edge_info_t ei = info.m_edge_map[e];
        edges_elems.push_back(get<1>(ei));
        edges_elems.push_back(get<2>(ei));
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
        h5h_write_dset(grp, "vertices", comm.m_vertices);
        h5h_write_dset(grp, "edges", comm.m_edges);
        h5h_write_dset(grp, "coherency", comm.m_coherency);
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
    if (nproc>1) mesh.partition_metis(nproc);

    string basename = fmesh;
    if (basename.find(".h5")  != string::npos) basename = basename.erase(basename.find(".h5"),  3);
    if (basename.find(".unv") != string::npos) basename = basename.erase(basename.find(".unv"), 4);
    mesh.write_proc_field(basename+".h5");

    for(int k=0;k<nproc;++k) {
        snprintf(fname, 1024,"%s.%04d.h5", argv[3], k);
        mesh.write_proc_file(fname, k);
    }
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
