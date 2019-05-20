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
#include "meshbase.h"
#include "aabb.h"
#include "sem_input.h" // sem_config_t

class Mesh3D
{
public:
    // methods
    Mesh3D():m_xadj(0L), m_adjncy(0L), debug(false), has_mrrs(false) {
	m_elems_offs.push_back(0);
    }

    void generate_output(int nprocs, const sem_config_t* config = NULL);
    void build_sf_interface();
    void compute_pml_free_surface();

    size_t n_nodes()     const { return m_elems.size(); }
    size_t n_vertices()  const { return m_xco.size(); }
    size_t n_elems()     const { return m_mat.size(); }
    size_t n_parts()     const { return n_procs; }
    size_t n_materials() const { return m_materials.size(); }
    size_t n_surfaces(const std::string c) const {
        size_t cpt=0;
        for (size_t i=0; i< m_surf_matname.size(); i++) {
            if (m_surf_matname[i].compare(0,c.length(),c)==0) cpt++;
        }
        return cpt;
    }
    size_t n_surface() const { return m_surf_matname.size();}

    void set_control_nodes(int n) { n_ctl_nodes = n; }

    int nodes_per_elem() const { return n_ctl_nodes; }

    index_t add_node(double x, double y, double z);
    index_t add_elem(int mat_idx, const Elem& el);
    index_t add_elem_mrrs(int mat_idx, int pos_mrrs, const Elem& el);

    int read_materials(const std::string& fname);
    int read_materials_v2(const std::string& fname);
    void write_materials(const std::string& fname);
    void write_materials_v2(const std::string& fname);
    void define_associated_materials();
    void read_mesh_file(const std::string& fname);
    void findelem(const std::vector<int>& elems, std::vector<int>& element, int & elnum) const;
    //size_t getData_line(char **buffer, size_t linesize, FILE* f);
    void partition_mesh(index_t n_parts);
    void dump_connectivity(const char* fname);

    index_t elem_part(index_t iel) const { return m_procs[iel]; }

    bool get_common_face(index_t e0, index_t e1, PFace& fc);

    void get_elem_nodes(index_t el, index_t nodes[8]) {
        index_t off = m_elems_offs[el];
        for(unsigned k=0;k<8;++k) nodes[k] = m_elems[off+k];
    }
    /// Returns domain number for an element:
    /// for now domain number == domain type:
    /// 1: fluid pml
    /// 2: solid pml
    /// 3: fluid
    /// 4: solid
    int get_elem_domain(index_t el) const {
        const Material& mat = m_materials[m_mat[el]];
        //printf("%d -> %c (%d/%d)\n", el, mat.ctype, m_mat[el], int(m_materials.size()));
        return mat.domain();
    }
    Surface* get_surface(const std::string& surfname) {
        Surface* sfp = m_surfaces[surfname];
        if (sfp!=NULL) return sfp;
        sfp = new Surface(surfname);
        m_surfaces[surfname] = sfp;
        return sfp;
    }

    void get_neighbour_elements(int nn, const index_t* n, std::set<index_t>& elemset) const;

public:
    // attributes
    index_t n_procs;
    int n_points;
    int n_neu;
    int n_PW;
    int ngrps;
    int n_ctl_nodes; ///< Number of control nodes per element (8 or 27)

    index_t *m_xadj, *m_adjncy;
    bool debug;
    bool has_mrrs;

    std::map<int,AABB> m_bbox;
    std::vector<double> m_xco,m_yco,m_zco;  ///< Coordinates of the nodes
    std::vector<index_t> m_elems; ///< size=(8|27)*n_elems ; describe each node of every elements
    //std::vector<int> m_Quad;
    std::vector<index_t> m_elems_offs; ///< size=n_elems+1; offset of node idx into elems
    std::vector<int> m_mat;  ///< size=n_elems; material index for element
    std::vector<int> m_mrrs;  ///< size=n_elems; mirror position for element
    std::vector<int> m_nelems_per_proc; // ?? number of elements for each procs
    std::vector<Material> m_materials;
    std::vector<unsigned int> m_vertex_domains;
    VertexElemMap  m_vertex_to_elem;
    std::vector<std::string> m_surf_matname;

    surf_info_map_t surfelem; // hexa8_num,<hexa8_id>,hexa8_tag,Quad4_tag

    void build_vertex_to_elem_map();
    void save_bbox();
    // A map of surfaces, indexed by names
    std::map<std::string,Surface*> m_surfaces;
protected:
    std::vector<index_t> m_procs; ///< size=n_elems; elem->proc association

protected:

    void read_mesh_hexa8(hid_t fid);
    void read_mesh_hexa27(hid_t fid);
    void read_mesh_Quad8(hid_t fid);
};


#endif

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
