/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// mesh.h : Gestion maillage format SEM
#ifndef _MESHPART_H_
#define _MESHPART_H_

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

/// Faces are stored as v1 v2 v3 v4, with v1 < v2 < v4
/// The nodes are consecutive, the direction is determined by the ordering of
/// the vertex number of the second node ie
/// if we create a face with (a,b,c,d) such that node numbers are c<b<a<d
/// then we store it as c,b,a,d so that with a global numbering of nodes
/// there is only one way to define the face and give it an orientation.
/// The domain to which the face belongs is added as a fifth component
/// because we want to duplicate the faces across domains (but keep them
/// oriented in the same manner)
struct PFace {
    PFace( int v[4], int dom ) {
        n[4] = dom;
        set_face(v);
    }
    PFace(const PFace& fc) { for(int k=0;k<5;++k) n[k]=fc.n[k]; }

    void set_face(int v[4]) {
        int l=0;
        int k;
        for(k=1;k<4;++k) {
            if (v[k]<v[l]) l=k;
        }
        int prev = (l+3)%4;
        int next = (l+1)%4;
        int step;
        if (v[prev]<v[next]) step=3;
        else step=1;
        for(k=0;k<4;++k) {
            n[k] = v[l];
            l = (l+step)%4;
        }
    }
    bool operator<(const PFace& fc) const {
	for(int i=0;i<5;++i) {
	    if (n[i] < fc.n[i]) return true;
	    if (n[i] > fc.n[i]) return false;
	}
	return false;
    }

    int n[5];
};

struct PEdge {
    PEdge( int v0, int v1, int dom ) {
        n[2] = dom;
        set_edge(v0, v1);
    }
    PEdge(const PEdge& ed) { for(int k=0;k<3;++k) n[k]=ed.n[k]; }

    void set_edge(int v0, int v1) {
        if (v0<v1) {
            n[0] = v0;
            n[1] = v1;
        } else {
            n[0] = v1;
            n[1] = v0;
        }
    }
    bool operator<(const PEdge& ed) const {
	for(int i=0;i<3;++i) {
	    if (n[i] < ed.n[i]) return true;
	    if (n[i] > ed.n[i]) return false;
	}
	return false;
    }

    int n[3];
};

typedef std::pair<int,int> PVertex; // pair(ID,Domain)

typedef std::map<PFace,int>  face_map_t;
typedef std::map<PEdge,int>  edge_map_t;
typedef std::map<PVertex,int> vertex_map_t;

/** Manages the list of communications with another processor */
class MeshPartComm
{
public:
    MeshPartComm() {}
    int m_dest; /// Destination processor

    /// The indices stored are in local numbering
    /// Before output, the local number must be reordered
    /// according to a global numbering scheme so that
    /// each communicating processor sees the items in the
    /// same order
    ///
    /// Note: there will be some duplication since if a face
    /// is shared, so are its edge and vertices.
    /// SEM can optimize this by not exchanging already registered
    /// gll
    std::vector<int> m_vertices;
    std::vector<int> m_edges;
    std::vector<int> m_faces;
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
    int add_facet(int n[4], int dom);
    int add_edge(int v0, int v1, int dom);
    int add_vertex(int v0, int dom);
    int add_node(int v0);

    int n_nodes() const    { return m_nodes_to_id.size(); }
    int n_elems() const    { return m_elems.size(); }
    int n_faces() const    { return m_face_to_id.size(); }
    int n_edges() const    { return m_edge_to_id.size(); }
    int n_vertices() const { return m_vertex_to_id.size(); }


    void get_local_nodes(std::vector<double>& nodes) const;
    void get_local_elements(std::vector<int>& elems) const;
    void get_local_materials(std::vector<int>& mats) const;
    void get_local_faces(std::vector<int>& mats) const;
    void get_local_edges(std::vector<int>& mats) const;
protected:
    int m_proc;
    const Mesh3D& m_mesh;

    std::vector<int> m_elems; // Local elements
    std::vector<int> m_elems_faces; // local face number of each local element
    std::vector<int> m_elems_edges; // local edge number of each local element
    std::vector<int> m_elems_vertices; // local vertex number of each local element
    face_map_t m_face_to_id;
    edge_map_t m_edge_to_id;
    vertex_map_t m_vertex_to_id;

    std::map<int,MeshPartComm> m_comm;
    std::map<int,int> m_nodes_to_id;

    void handle_local_element(int el);
    void handle_neighbour_element(int el);
};

#endif

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
