/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include <cstdio>
#include "mesh.h"
#include "metis.h"
#include <map>
#include <cstdlib>

using std::map;
using std::multimap;
using std::vector;
using std::pair;

// =====================================================





// =====================================================
int Mesh3D::add_node(double x, double y, double z)
{
    m_xco.push_back( x );
    m_yco.push_back( y );
    m_zco.push_back( z );
    return m_xco.size()-1;
}

int Mesh3D::add_elem(int mat_idx, const HexElem& el)
{
    // Builds elem<->vertex graph
    for(int i=0;i<8;++i) {
	m_elems.push_back(el.v[i]);
    }
    m_elems_offs.push_back(m_elems.size());
    m_mat.push_back( mat_idx );
    return m_elems_offs.size()-1;
}


void Mesh3D::partition_mesh(int n_parts)
{
    int ne = n_elems();
    int nn = n_vertices();
    int ncommon = 1;
    int numflags = 0;
    int ncon=1;
    int *vwgt=0L;
    int *vsize=0L;
    int *adjwgt=0L;
    float *tpwgts=0L;
    float *ubvec=0L;
    int *options=0L;
    int edgecut;

    n_procs = n_parts;
    m_procs.resize(ne);

    METIS_MeshToDual(&ne, &nn, &m_elems_offs[0], &m_elems[0],
		     &ncommon, &numflags, &m_xadj, &m_adjncy);

    if (n_parts>1) {
        METIS_PartGraphKway(&ne, &ncon, m_xadj, m_adjncy,
                            vwgt, vsize, adjwgt, &n_procs, tpwgts, ubvec,
                            options, &edgecut, &m_procs[0]);
    } else {
        for(int k=0;k<ne;++k) m_procs[k]=0;
    }
}





int Mesh3D::read_materials(const std::string& str)
{
    FILE* f = fopen(str.c_str(), "r");
    int nmats;
    char type;
    double vs, vp, rho;
    int ngllx, nglly, ngllz;
    double dt, Qp, Qmu;

    fscanf(f, "%d", &nmats);
    for(int k=0;k<nmats;++k) {
        fscanf(f, "%c %lf %lf %lf %d %d %d %lf %lf %lf",
               &type, &vs, &vp, &rho, &ngllx, &nglly, &ngllz,
               &dt, &Qp, &Qmu);
        m_materials.push_back(Material(type, vp, vs, rho, Qp, Qmu, ngllx, nglly, ngllz));
    }
    return nmats;
}

int Mesh3D::add_material()
{
    m_materials.push_back(Material());
    return m_materials.size()-1;
}

void Mesh3D::read_mesh_file(const std::string& fname)
{
    hid_t file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    h5h_read_dset_Nx3(file_id, "/Nodes", m_xco, m_yco, m_zco);

    if (H5Lexists(file_id, "/Sem3D/Hexa8", H5P_DEFAULT)) {
        read_mesh_hexa8(file_id);
    } else if (H5Lexists(file_id, "/Sem3D/Hexa27", H5P_DEFAULT)) {
        read_mesh_hexa27(file_id);
    }
    else {
        printf("ERR: only Quad4 and Quad8 are supported\n");
        exit(1);
    }

    h5h_read_dset(file_id, "/Sem3D/Mat", m_mat);
}

void Mesh3D::read_mesh_hexa8(hid_t file_id)
{
    int nel, nnodes;
    h5h_read_dset_2d(file_id, "/Sem3D/Hexa8", nel, nnodes, m_elems);
    n_ctl_nodes = 8;
    if (nnodes!=8) {
        printf("Error: dataset /Sem/Hexa8 is not of size NEL*8\n");
        exit(1);
    }
    for(int k=0;k<nel;++k) {
        m_elems_offs.push_back(8*(k+1));
    }
}

void Mesh3D::read_mesh_hexa27(hid_t file_id)
{
}

void Mesh3D::build_vertex_to_elem_map()
{
    int nel = n_elems();
    m_vertex_to_elem.init(nel);
    for(int i=0;i<nel;++i) {
        for(int k=m_elems_offs[i];k<m_elems_offs[i+1];++k) {
            m_vertex_to_elem.add_link(m_elems[k], i);
        }
    }
}


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
