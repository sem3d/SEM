/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// mesh.h : Gestion maillage format SEM
#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <string>
#include <cstdio>
#include <map>
#include "material.h"
#include "h5helper.h"
#include "vertex_elem_map.h"

template <int N>
struct Elem {
    Elem() {}
    Elem(const Elem& el) { for(int i=0;i<N;++i) { v[i] = el.v[i];} }
    Elem& operator=(const Elem& el) { for(int i=0;i<N;++i) { v[i] = el.v[i]; } return *this; }
    bool operator==(const Elem& el) {
	for(int i=0;i<N;++i) { if (v[i] != el.v[i]) return false; }
	return true;
    }
    bool operator<(const Elem& el) {
	for(int i=0;i<N;++i) {
	    if (v[i] < el.v[i]) return true;
	    if (v[i] > el.v[i]) return false;
	}
	return false;
    }
    int v[N];
};

struct HexElem : public Elem<8>
{
    HexElem() {}
    HexElem(int a, int b, int c, int d,
	    int e, int f, int g, int h) {
	v[0] = a;
	v[1] = b;
	v[2] = c;
	v[3] = d;
	v[4] = e;
	v[5] = f;
	v[6] = g;
	v[7] = h;
    }
};

struct QuadElem : public Elem<4>
{
    QuadElem() {}
    QuadElem(int a, int b, int c, int d){
	v[0] = a;
	v[1] = b;
	v[2] = c;
	v[3] = d;
    }
};

struct FaceDesc {
    int v[4];  /// Local vertex index
    int e[4];  /// Local edge index

    void show_face() {
	printf("%d.%d.%d.%d\n", v[0], v[1], v[2], v[3]);
    }
};

class Mesh3D
{
public:
    // methods
    Mesh3D():m_xadj(0L), m_adjncy(0L) {
	m_elems_offs.push_back(0);
    }

    int n_nodes()     const { return m_elems.size(); }
    int n_vertices()  const { return m_xco.size(); }
    int n_elems()     const { return m_mat.size(); }
    int n_parts()     const { return n_procs; }
    int n_materials() const { return m_materials.size(); }

    int nodes_per_elem() const { return 8; }

    int add_node(double x, double y, double z);
    int add_elem(int mat_idx, const HexElem& el);

    int read_materials(const std::string& fname);
    void read_mesh_file(const std::string& fname);


    void partition_mesh(int n_procs);

    int elem_part(int iel) const { return m_procs[iel]; }

    int n_shared_faces() const { return -1; }
    int n_shared_edges() const { return -1; }
    int n_shared_vertices() const { return -1; }

    void get_elem_nodes(int el, int nodes[8]) {
        int off = m_elems_offs[el];
        for(int k=0;k<7;++k) nodes[k] = m_elems[off+k];
    }
    /// Returns domain number for an element:
    /// for now domain number == domain type:
    /// 1: fluid pml
    /// 2: solid pml
    /// 3: fluid
    /// 4: solid
    int get_elem_domain(int el) const {
        const Material& mat = m_materials[m_mat[el]];
        //printf("%d -> %c (%d/%d)\n", el, mat.ctype, m_mat[el], int(m_materials.size()));
        return mat.domain();
    }
public:
    // attributes
    int n_procs;
    int n_points;
    int n_neu;
    int n_PW;
    int n_ctl_nodes; ///< Number of control nodes per element (8 or 27)

    int *m_xadj, *m_adjncy;

    std::vector<double> m_xco,m_yco,m_zco;  ///< Coordinates of the nodes
    std::vector<int> m_elems; ///< size=8*n_elems ; describe each node of every elements
    std::vector<int> m_elems_offs; ///< size=n_elems+1; offset of node idx into elems
    std::vector<int> m_mat;  ///< size=n_elems; material index for element
    std::vector<int> m_nelems_per_proc; // ?? number of elements for each procs
    std::vector<Material> m_materials;
    VertexElemMap  m_vertex_to_elem;
    void build_vertex_to_elem_map();
protected:
    std::vector<int> m_procs; ///< size=n_elems; elem->proc association

protected:

    void read_mesh_hexa8(hid_t fid);
    void read_mesh_hexa27(hid_t fid);
};


#endif

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
