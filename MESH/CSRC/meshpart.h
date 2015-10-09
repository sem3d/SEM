/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// mesh.h : Gestion maillage format SEM
#ifndef _MESHPART_H_
#define _MESHPART_H_

#include "meshbase.h"
#include "mesh.h"
#include <vector>
#include <map>


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

    void get_faces(std::vector<int>& tmp) const {
        tmp.clear();
        for(face_map_t::const_iterator it=m_faces.begin();it!=m_faces.end();++it) {
            tmp.push_back(it->second);
        }
    }
    void get_edges(std::vector<int>& tmp) const {
        tmp.clear();
        for(edge_map_t::const_iterator it=m_edges.begin();it!=m_edges.end();++it) {
            tmp.push_back(it->second);
        }
    }
    void get_vertices(std::vector<int>& tmp) const {
        tmp.clear();
        for(vertex_map_t::const_iterator it=m_vertices.begin();it!=m_vertices.end();++it) {
            tmp.push_back(it->second);
        }
    }
};

/** Manages a part of the mesh on a specific processor */
class Mesh3DPart {
public:
    Mesh3DPart(const Mesh3D& mesh, int proc):
        m_mesh(mesh),
        m_proc(proc) {}


    void compute_part();
    void output_mesh_part();
    void output_mesh_part_xmf();

    // manages local facets
    int add_facet(const PFace& facet);
    int add_edge(const PEdge& edge);
    int add_vertex(const PVertex& vertex);
    int add_facet(int n[4], int dom) { PFace facet(n, dom); return add_facet(facet); }
    int add_edge(int v0, int v1, int dom) { PEdge edge(v0, v1, dom); return add_edge(edge); }
    int add_vertex(int v0, int dom) { PVertex vertex(v0, dom); return add_vertex(vertex); }
    int add_node(int v0);

    int n_nodes() const    { return m_nodes_to_id.size(); }
    int n_elems() const    { return m_elems.size(); }
    int n_faces() const    { return m_face_to_id.size(); }
    int n_edges() const    { return m_edge_to_id.size(); }
    int n_vertices() const { return m_vertex_to_id.size(); }


    void get_local_nodes(std::vector<double>& nodes) const;
    void get_local_elements(std::vector<int>& elems) const;
    void get_local_materials(std::vector<int>& mats, std::vector<int>& doms) const;
    void get_local_faces(std::vector<int>& faces, std::vector<int>& doms) const;
    void get_local_edges(std::vector<int>& edges, std::vector<int>& doms) const;
    void get_face_coupling(int d0, int d1, std::vector<int>& cpl, std::vector<int>& orient) const;
    void get_edge_coupling(int d0, int d1, std::vector<int>& cpl) const;
    void get_vertex_coupling(int d0, int d1, std::vector<int>& cpl) const;

    // Generic function to write and compute coupling interfaces
    void write_coupling_interface(hid_t fid, const char* pfx, int d0, int d1);
    void output_int_scalar(FILE* f, int indent, const char* aname,
                           const char* atype, int n0, const char* field);
    void reorder_comm(MeshPartComm& comm);
protected:
    const Mesh3D& m_mesh;
    int m_proc;

    std::vector<int> m_elems; // Local elements
    std::vector<int> m_elems_faces; // local face number of each local element
    std::vector<int> m_elems_edges; // local edge number of each local element
    std::vector<int> m_elems_vertices; // local vertex number of each local element
    std::vector<Surface*> m_surfaces;
    face_map_t m_face_to_id;
    edge_map_t m_edge_to_id;
    vertex_map_t m_vertex_to_id;

    std::map<int,MeshPartComm> m_comm;
    std::map<int,int> m_nodes_to_id;

    void handle_local_element(int el);
    void handle_neighbour_element(int el);
    void handle_surface(const Surface* surf);
    void output_mesh_attributes(hid_t fid);
    void output_local_mesh(hid_t fid);
    void output_surface(hid_t fid, const Surface* surf);
    void write_surface_dom(hid_t gid, const Surface* surf, const char* pfx, int dom);
};

#endif

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
