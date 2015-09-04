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


struct edge_t {
    edge_t() {}
    edge_t(const edge_t& e):v0(e.v0), v1(e.v1),n(e.n),elems(e.elems) {}
    edge_t(int _v0, int _v1, int _n, int _el):v0(_v0),v1(_v1),n(_n)
	{ elems.push_back(_el); }

    bool operator<(const edge_t& e) const {
	if (v0<e.v0) return true;
	if (v0==e.v0) return v1<e.v1;
	return false;
    }
    bool operator==(const edge_t& e) const {
	return v0==e.v0 && v1==e.v1;
    }

    int orient(int ve0, int ve1) {
	if ((ve0==v0)&&(ve1==v1)) return 0;
	if ((ve0==v0)&&(ve1==v1)) return 1;
	return -1; // Error
    }
    int v0, v1;
    int n; // Edge number
    std::vector<int> elems; // elements whose this edge belongs to
};

struct ordered_edge_t
{
    ordered_edge_t(const ordered_edge_t& e):v0(e.v0), v1(e.v1) {}
    ordered_edge_t(int _v0, int _v1):v0(_v0),v1(_v1) {
	if (_v1<_v0) {
	    v0 = _v1;
	    v1 = _v0;
	}
    }
    bool operator<(const ordered_edge_t& e) const {
	if (v0<e.v0) return true;
	if (v0==e.v0) return v1<e.v1;
	return false;
    }
    bool operator==(const ordered_edge_t& e) const {
	return v0==e.v0 && v1==e.v1;
    }
    int v0, v1;
};

/// Contains information about structures (face, edge, vertice) shared between processor
/// The structure for p0->p1 should match p1->p0 ie :
/// foreach i p0->local_to_global(p0->vertices[i])=p1->local_to_global(p1->vertices[i])
struct shared_part_t
{
    std::vector<int> vertices;
    std::vector<int> faces;
    std::vector<int> faces_orient;
    std::vector<int> edges;
    std::vector<int> edges_orient;
};

/**
  This class holds the parts of a mesh
  local to one processor
*/
class MeshPart
{
public:
    int n_elems()    { return n_elems_per_proc; }
    int n_faces()    { return m_n_faces; }
    int n_edges()    { return m_n_edges; }
    /// Total number of distinct vertices (always 8 for one hex)
    int n_vertices() { return n_l2g.size(); }
    /// Number of control nodes (8 or 27 for one hex)
    int n_nodes()    { return n_l2g.size(); }

    int global_node_idx(int k) { return n_l2g[k]; }
public:
    int part;
    int n_elems_per_proc;
    int m_n_faces;
    int m_n_edges;
    std::vector<int> n_l2g;
    std::vector<int> e_l2g;
    std::vector<int> n_g2l;
    std::vector<int> e_g2l;

    std::vector<int> faces;
    std::vector<int> faces_orient;

    std::vector<int> edges;
    std::vector<int> edges_orient;

    std::vector<int> vertices;

};

struct FaceDesc {
    int v[4];  /// Local vertex index
    int e[4];  /// Local edge index

    void show_face() {
	printf("%d.%d.%d.%d\n", v[0], v[1], v[2], v[3]);
    }
};

// A compact graph structure.
struct graph_t
{
    graph_t() { offsets.push_back(0); }
    int n_elems() { return offsets.size()-1; }
    std::vector<int> links;
    std::vector<int> offsets;
};

class Mesh3D
{
public:
    // methods
    Mesh3D():m_xadj(0L), m_adjncy(0L) {
	m_elems_offs.push_back(0);
    }

    int n_nodes() { return m_elems.size(); }
    int n_vertices() { return m_xco.size(); }
    int n_elems() { return m_mat.size(); }
    int n_parts() { return n_procs; }
    int n_materials() { return m_materials.size(); }

    int nodes_per_elem() { return 8; }

    int add_node(double x, double y, double z);
    int add_elem(int mat_idx, const HexElem& el);
    int add_material(); // Aucune propriete pour l'instant

    int read_materials(const std::string& fname);
    void read_mesh_file(const std::string& fname);


    void partition_mesh(int n_procs);

    int elem_part(int iel) { return m_procs[iel]; }

    void compute_local_part(int part, MeshPart& loc);

    /// Accessors/reordering for local parts
    void get_local_material(MeshPart& loc, std::vector<int>& tmp);
    void get_local_nodes(MeshPart& loc, std::vector<double>& tmp);
    void get_local_faces(MeshPart& loc, std::vector<int>& tmp);
    void get_local_elements(MeshPart& loc, std::vector<int>& tmp);

    int n_shared_faces() const { return -1; }
    int n_shared_edges() const { return -1; }
    int n_shared_vertices() const { return -1; }

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
    std::map< ordered_edge_t, edge_t > m_edge_map; // Maps edges' vertices to edge definition
    std::multimap<int,int> m_elem_proc_map;
    std::map<std::pair<int,int>, shared_part_t> m_shared;
protected:
    std::vector<int> m_procs; ///< size=n_elems; elem->proc association

protected:
    /// Compute a list of elements sharing a border with another processor
    void compute_comm_elements();
    void compute_shared_objects_for_part(int part);
    void compute_nodes_indexes(MeshPart& loc);
    void compute_local_connectivity(MeshPart& loc);
    /// Adds a local face identified by it's elements and face number. Faces are shared between elements.
    void add_local_face(MeshPart& loc, int iel, int gel, int nf);
    void add_local_edge(MeshPart& loc, int iel, int gel, int ne);

    /// returns false if global elem el0 doesn't share it's local face nf0 with global el1
    /// if true returns the local face number of el1 and their relative orientation
    bool shared_face(int el0, int nf0, int el1, FaceDesc& other);

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
