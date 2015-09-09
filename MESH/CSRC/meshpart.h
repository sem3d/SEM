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


/// Faces are stored as v1 v2 v3 v4, with v1 < v2 < v4
/// The nodes are consecutive, the direction is determined by the ordering of
/// the vertex number of the second node ie
/// if we create a face with (a,b,c,d) such that node numbers are c<b<a<d
/// then we store it as c,b,a,d so that with a global numbering of nodes
/// there is only one way to define the face and give it an orientation.
/// The domain to which the face belongs is added as a fifth component
/// because we want to duplicate the faces across domains (but keep them
/// oriented in the same manner)
struct Face {
    Face( int v[4], int dom ) {
        n[4] = dom;
        set_face(v);
    }
    Face(const Face& fc) { for(int k=0;k<5;++k) n[k]=fc.n[k]; }

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
    bool operator<(const Face& fc) const {
	for(int i=0;i<5;++i) {
	    if (n[i] < fc.n[i]) return true;
	    if (n[i] > fc.n[i]) return false;
	}
	return false;
    }

protected:
    int n[5];
};


typedef std::map<Face,int>  face_map_t;

/** Manages a part of the mesh on a specific processor */
class Mesh3DPart {
public:
    Mesh3DPart(const Mesh3D& mesh, int proc):
        m_mesh(mesh),
        m_proc(proc) {}


    void compute_part();
    void output_mesh_part();
    void output_mesh_part_xmf();

    int add_facet(int n[4], int dom);
protected:
    int m_proc;
    const Mesh3D& m_mesh;

    face_map_t m_face_to_face_id;

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
