/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// meshbase.h : Elements de base d'un maillage hexa (face,edge,vertex)
#ifndef _MESHBASE_H_
#define _MESHBASE_H_

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <cstdlib>
#include <cassert>
#include "sem_types.h"

// type for global indexes
constexpr index_t invalid_index = (index_t)(-1);
// type for local indexes (to each cpu partition inmem)
typedef int32_t  loc_index_t;



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
    PFace() {}
    PFace( index_t v[4], int _dom) {
        set_domain(_dom);
        set_face(v);
    }
    void set_domain(int _dom) {
        dom = _dom;
    }
    PFace(const PFace& fc):dom(fc.dom),orient(fc.orient)
        { for(int k=0;k<4;++k) n[k]=fc.n[k]; }

    void set_face(index_t v[4]) {
        int l=0;
        int k;
        for(k=1;k<4;++k) {
            if (v[k]<v[l]) l=k;
        }
        int prev = (l+3)%4;
        int next = (l+1)%4;
        int step;
        if (v[prev]<v[next]) {
            step=3;
            orient = -1;
        } else {
            step=1;
            orient = 1;
        }
        for(k=0;k<4;++k) {
            n[k] = v[l];
            l = (l+step)%4;
        }
    }
    bool operator<(const PFace& fc) const {
	for(int i=0;i<4;++i) {
	    if (n[i] < fc.n[i]) return true;
	    if (n[i] > fc.n[i]) return false;
	}
        if (dom<fc.dom) return true;
	return false;
    }
    bool operator==(const PFace& fc) const {
	for(int i=0;i<4;++i) {
	    if (n[i] != fc.n[i]) return false;
	}
        if (dom!=fc.dom) return false;
	return true;
    }
    /// Same as operator= but ignore domain
    bool eq_geom(const PFace& fc) const {
	for(int i=0;i<4;++i) {
	    if (n[i] != fc.n[i]) return false;
	}
	return true;
    }
    int domain() const { return dom; }
    index_t n[4];
    int dom;
    // !! This will be only useful for faces at domain interfaces, otherwise
    // there is no way to tell which of the two elements sharing the face appears first
    int orient; // 1 if points inside original element -1 otherwise
};

struct PEdge {
    PEdge() {}
    PEdge( index_t v0, index_t v1, int dom ) {
        n[2] = dom;
        set_edge(v0, v1);
    }
    PEdge(const PEdge& ed) {
        for(int k=0;k<3;++k) n[k]=ed.n[k];
    }

    void set_edge(index_t v0, index_t v1) {
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
    bool operator==(const PEdge& ed) const {
	for(int i=0;i<3;++i) {
	    if (n[i] != ed.n[i]) return false;
	}
	return true;
    }
    /// Same as operator= but ignore domain
    bool eq_geom(const PEdge& ed) const {
	for(int i=0;i<2;++i) {
	    if (n[i] != ed.n[i]) return false;
	}
	return true;
    }
    int domain() const { return n[2]; }
    index_t n[3];
};


struct PVertex {
    PVertex() {}
    PVertex( index_t v0, int dom ) {
        n[0] = v0;
        n[1] = dom;
    }
    PVertex(const PVertex& vx) {
        for(int k=0;k<2;++k) n[k]=vx.n[k];
    }

    bool operator<(const PVertex& vx) const {
	for(int i=0;i<2;++i) {
	    if (n[i] < vx.n[i]) return true;
	    if (n[i] > vx.n[i]) return false;
	}
	return false;
    }
    bool operator==(const PVertex& vx) const {
	for(int i=0;i<2;++i) {
	    if (n[i] != vx.n[i]) return false;
	}
	return true;
    }
    /// Same as operator= but ignore domain
    bool eq_geom(const PVertex& vx) const {
        if (n[0] != vx.n[0]) return false;
	return true;
    }
    int domain() const { return n[1]; }
    index_t n[2];
};

struct surf_info_t {
    std::vector<index_t> xxnodes;
    int  dom;
    int  mat;
};

typedef std::map<PFace,index_t>  face_map_t;
typedef std::map<PEdge,index_t>  edge_map_t;
typedef std::map<PVertex,index_t> vertex_map_t;
typedef std::map<index_t,index_t> node_id_map_t;
typedef std::map<index_t, surf_info_t  > surf_info_map_t;

class Surface {
public:
    Surface(const std::string& name):m_name(name) {}
    Surface(const Surface& surf):m_faces(surf.m_faces) {}

    bool empty() const { return m_faces.size()==0; }
    void add_face(const PFace& fc, int data) {
        m_faces[fc] = data;
    }
    void add_edge(const PEdge& ed, int data) {
        m_edges[ed] = data;
    }
    void add_vertex(const PVertex& vx, int data) {
        m_vertices[vx] = data;
    }
    void get_faces_data( int dom, std::vector<int>& data,
                         std::vector<int>& orient,
                         std::vector<int>& matdom) const {
        data.clear();
        orient.clear();
        matdom.clear();

        for(face_map_t::const_iterator it=m_faces.begin();it!=m_faces.end();++it) {
            if (it->first.domain()!=dom) continue;
            data.push_back(it->second);
            orient.push_back(it->first.orient);
            matdom.push_back(it->first.domain());
        }
    }
    void get_edges_data(int dom, std::vector<int>& data, std::vector<int>& orient, std::vector<int>& matdom) const {
        data.clear();
        matdom.clear();
        orient.clear();
        for(edge_map_t::const_iterator it=m_edges.begin();it!=m_edges.end();++it) {
            if (it->first.domain()!=dom) continue;
            data.push_back(it->second);
            matdom.push_back(it->first.domain());
        }
    }
    void get_vertices_data(int dom, std::vector<int>& data, std::vector<int>& matdom) const {
        data.clear();
        matdom.clear();
        for(vertex_map_t::const_iterator it=m_vertices.begin();it!=m_vertices.end();++it) {
            if (it->first.domain()!=dom) continue;
            data.push_back(it->second);
            matdom.push_back(it->first.domain());
        }
    }
    std::string name() const { return m_name; }
    
    std::string m_name;
    face_map_t m_faces;
    edge_map_t m_edges;
    vertex_map_t m_vertices;
};


struct Elem {
    Elem(int _N):N(_N) { v=(int*)malloc(N*sizeof(int)); }
    ~Elem() { free(v); }
    Elem(const Elem& el):N(el.N) { v=(int*)malloc(N*sizeof(int)); for(int i=0;i<N;++i) { v[i] = el.v[i];} }
    Elem& operator=(const Elem& el) { assert(N==el.N); for(int i=0;i<N;++i) { v[i] = el.v[i]; } return *this; }
    bool operator==(const Elem& el) {
        assert(N==el.N);
	for(int i=0;i<N;++i) { if (v[i] != el.v[i]) return false; }
	return true;
    }
    bool operator<(const Elem& el) {
        assert(N==el.N);
	for(int i=0;i<N;++i) {
	    if (v[i] < el.v[i]) return true;
	    if (v[i] > el.v[i]) return false;
	}
	return false;
    }
    int N;
    int *v;
};

struct HexElem : public Elem
{
    HexElem(int nn=8):Elem(nn) {}
};


struct FaceDesc {
    int v[4];  /// Local vertex index
    int e[4];  /// Local edge index

    void show_face() {
	printf("%d.%d.%d.%d\n", v[0], v[1], v[2], v[3]);
    }
};

/// Reference numbering for faces
/// Faces are defined so that the normal given using the right hand rule
/// points inside
/// This definition needs to match indexation.f90
static const FaceDesc RefFace[6] = {
    { {0,1,2,3}, {0, 1, 2, 3} },
    { {0,4,5,1}, {0, 4, 5, 6} },
    { {1,5,6,2}, {1, 7, 8, 4} },
    { {3,2,6,7}, {2, 7, 9,10} },
    { {0,3,7,4}, {3,10,11, 6} },
    { {4,7,6,5}, {5, 8, 9,11} },
};

/// Reference numbering of edges
static const int RefEdge[12][2] = {
    { 0, 1},
    { 1, 2},
    { 3, 2},
    { 0, 3},
    { 4, 5},
    { 5, 6},
    { 7, 6},
    { 4, 7},
    { 0, 4},
    { 1, 5},
    { 2, 6},
    { 3, 7}
};

#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
