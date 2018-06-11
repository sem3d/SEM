/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// mesh.h : Gestion maillage format SEM
#ifndef _MESHPART_H_
#define _MESHPART_H_

#include "meshbase.h"
#include "mesh.h"
#include "sem_input.h" // sem_config_t
#include <vector>
#include <map>
#include <set>

template<typename T0, typename T1>
T1 get(const std::map<T0,T1>& m, const T0& key, const T1& def)
{
    typename std::map<T0,T1>::const_iterator it;
    it = m.find(key);
    if (it!=m.end()) return it->second;
    return def;
}


/** Manages the list of communications with another processor */
class MeshPartComm
{
public:
    MeshPartComm():m_dest(-1) {
    }
    int m_dest; /// Destination processor

    /// Note: there will be some duplication since if a face
    /// is shared, so are its edge and vertices.
    /// SEM can optimize this by not exchanging already registered
    /// gll
    ///
    /// we use a map to store the shared face/edge/vertex
    /// since std::map maintains its keys in order,
    /// it ensures they are ordered the same way on both cpus
    vertex_map_t m_vertices;
    edge_map_t m_edges;
    face_map_t m_faces;

    size_t n_faces() const { return m_faces.size(); }
    size_t n_edges() const { return m_edges.size(); }
    size_t n_vertices() const { return m_vertices.size(); }

    void get_faces(std::vector<index_t>& defs,
                   std::vector<index_t>& ids, std::vector<int>& doms) const {
        defs.clear();
        ids.clear();
        doms.clear();
        for(face_map_t::const_iterator it=m_faces.begin();it!=m_faces.end();++it) {
            for(int k=0;k<4;++k) defs.push_back(it->first.n[k]);
            ids.push_back(it->second);
            doms.push_back(it->first.domain());
        }
    }
    void get_edges(std::vector<index_t>& defs,
                   std::vector<index_t>& ids, std::vector<int>& doms) const {
        defs.clear();
        ids.clear();
        doms.clear();
        for(edge_map_t::const_iterator it=m_edges.begin();it!=m_edges.end();++it) {
            for(int k=0;k<2;++k) defs.push_back(it->first.n[k]);
            ids.push_back(it->second);
            doms.push_back(it->first.domain());
        }
    }
    void get_vertices(std::vector<index_t>& defs,
                      std::vector<index_t>& ids,
                      std::vector<int>& doms) const {
        defs.clear();
        ids.clear();
        doms.clear();
        for(vertex_map_t::const_iterator it=m_vertices.begin();it!=m_vertices.end();++it) {
            defs.push_back(it->first.n[0]);
            ids.push_back(it->second);
            doms.push_back(it->first.domain());
        }
    }
};

/** Manages a part of the mesh on a specific processor */
class Mesh3DPart {
public:
    Mesh3DPart(const Mesh3D& mesh, int proc, const sem_config_t* config):
        m_mesh(mesh),
        m_proc(proc),
        m_cfg(config) {}


    void compute_part();
    void compute_face_communications();
    void compute_edge_communications();
    void compute_vertex_communications();
    void output_mesh_part();
    void output_mesh_part_xmf();
    /// Returns whether an element (global index) communicates with other processors
    bool is_border_element(int el);

    // manages local facets
    index_t add_facet(const PFace& facet, bool border);
    index_t add_edge(const PEdge& edge, bool border);
    index_t add_vertex(const PVertex& vertex, bool border);
    index_t add_node(index_t v0);

    int get_mat_(int el) const;
    size_t n_nodes() const    { return m_nodes_to_id.size(); }
    size_t n_elems() const    { return m_elems.size(); }
    size_t n_faces() const    { return m_face_to_id.size(); }
    size_t n_edges() const    { return m_edge_to_id.size(); }
    size_t n_vertices() const { return m_vertex_to_id.size(); }

    void define_Eventual_Boundary_Surface(int el);
    void get_local_nodes(std::vector<double>& nodes) const;
    void get_local_elements(std::vector<loc_index_t>& elems) const;
    void get_local_materials(std::vector<int>& mats, std::vector<int>& doms) const;
    void get_local_mirrors(std::vector<int>& pmrrs) const;
    void get_local_faces(std::vector<loc_index_t>& faces, std::vector<int>& doms) const;
    void get_local_edges(std::vector<loc_index_t>& edges, std::vector<int>& doms) const;

    void get_local_vertices_dom(std::vector<int>& doms) const;
    void get_face_coupling(int d0, int d1, std::vector<index_t>& cpl, std::vector<int>& orient) const;
    void get_edge_coupling(int d0, int d1, std::vector<index_t>& cpl) const;
    void get_vertex_coupling(int d0, int d1, std::vector<index_t>& cpl) const;

    // Generic function to write and compute coupling interfaces
    void write_coupling_interface(hid_t fid, const char* pfx, int d0, int d1);
    void output_int_scalar(FILE* f, int indent, const char* aname,
                           const char* atype, int n0, const char* field);
    void output_int_constant(FILE* f, int indent, const char* aname, const char* atype, int val);
    void reorder_comm(MeshPartComm& comm);
    void global_to_local_ids(std::vector<index_t>& ids) const;

protected:
    const Mesh3D& m_mesh;
    int m_proc;
    const sem_config_t* m_cfg;

    std::vector<index_t> m_elems; // Local elements
    std::vector<index_t> m_elems_faces; // local face number of each local element
    std::vector<index_t> m_elems_edges; // local edge number of each local element
    std::vector<index_t> m_elems_vertices; // local vertex number of each local element
    std::vector<Surface*> m_surfaces;
    face_map_t m_face_to_id;
    edge_map_t m_edge_to_id;
    vertex_map_t m_vertex_to_id;
    // Maintains a flag for whose element might be in contact with other processors
    std::vector<bool> m_face_border;
    std::vector<bool> m_edge_border;
    std::vector<bool> m_vertex_border;

    std::map<int,MeshPartComm> m_comm;
    node_id_map_t m_nodes_to_id;

    // Mirror
    std::vector<index_t> m_mirror_e;
    std::vector<index_t> m_mirror_ijk;
    std::vector<double> m_mirror_xyz;
    std::vector<double> m_mirror_w;
    std::vector<double> m_mirror_inside; // Inside <=> window [0., 1.]
    std::vector<double> m_mirror_outnormal;
    std::vector<double> m_gll;
    void compute_gll();
    void handle_mirror_implicit_surf(index_t el);
    void handle_mirror_surf(const Surface* smirror);
    void shape8_local2global(double const vco[3][8],
                             const double& xi, const double& eta, const double& zeta,
                             double& x, double& y, double& z) const;

    void handle_local_element(index_t el, bool is_border);
    void handle_neighbouring_face(index_t lnf, const PFace& fc, index_t el);
    void handle_neighbouring_edge(index_t lne, const PEdge& ed, index_t el);
    void handle_neighbouring_vertex(index_t lnv, const PVertex& vx, index_t el);
    void handle_surface(const Surface* surf);
    void output_mesh_attributes(hid_t fid);
    void output_local_mesh(hid_t fid);
    void output_surface(hid_t fid, const Surface* surf);
    void output_mirror() const;
    void output_comm(hid_t gid, const MeshPartComm& comm, int dest);

    void write_surface_dom(hid_t gid, const Surface* surf, const char* pfx, int dom);
    void get_neighbour_elements(int nn, const int* n, std::set<int>& elemset);
    void output_xmf_elements();
    void output_xmf_faces();
    void output_xmf_edges();
    void output_xmf_vertices();
    void output_xmf_mirror();
    void output_xmf_comms();
    void output_xmf_comms_faces();
    void output_xmf_comms_edges();
    void output_xmf_header(FILE* f);
    void output_xmf_footer(FILE* f);
};

#endif

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
