/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include <cstdio>
#include "mesh.h"
#include "metis.h"
#include <map>
#include <cstdlib>
#include <cstring>

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
    vector<int> vwgt;
    int *vsize=0L;
    int *adjwgt=0L;
    float *tpwgts=0L;
    float *ubvec=0L;
    int *options=0L;
    int edgecut;

    n_procs = n_parts;
    m_procs.resize(ne);
    m_xadj = 0L;
    m_adjncy = 0L;
    METIS_MeshToDual(&ne, &nn, &m_elems_offs[0], &m_elems[0],
		     &ncommon, &numflags, &m_xadj, &m_adjncy);

    //dump_connectivity("conn1.dat");
    // Tentative de reordonnancement des elements pour optimiser la reutilisation de cache
    // lors de la boucle sur les elements
//    vector<int> perm, iperm;
//    perm.resize(ne);
//    iperm.resize(ne);
//    METIS_NodeND(&ne, m_xadj, m_adjncy, 0L, 0L, &perm[0], &iperm[0]);
//    for(int k=0;k<m_xadj[ne];++k) {
//        m_adjncy[k] = perm[m_adjncy[k]];
//    }
//    dump_connectivity("conn2.dat");
    vwgt.resize(ne);
    // Define weights
    for(int k=0;k<ne;++k) {
        const Material& mat = m_materials[m_mat[k]];
        switch(mat.m_type) {
        case DM_SOLID:
            vwgt[k] = 3;
            break;
        case DM_FLUID:
            vwgt[k] = 1;
            break;
        case DM_SOLID_PML:
            vwgt[k] = 9;
            break;
        case DM_FLUID_PML:
            vwgt[k] = 3;
            break;
        default:
            vwgt[k] = 1;
        }
    }
    if (n_parts>1) {
        METIS_PartGraphKway(&ne, &ncon, m_xadj, m_adjncy,
                            &vwgt[0], vsize, adjwgt, &n_procs, tpwgts, ubvec,
                            options, &edgecut, &m_procs[0]);
    } else {
        for(int k=0;k<ne;++k) m_procs[k]=0;
    }
}


void Mesh3D::dump_connectivity(const char* fname)
{
    FILE* fmat = fopen(fname, "wb");
    int ne = n_elems();
    unsigned char* mat = (unsigned char*)malloc(ne*ne*sizeof(unsigned char));
    memset(mat, 0, ne*ne);
    for(int i=0;i<n_elems();++i) {
        for(int k=m_xadj[i];k<m_xadj[i+1];++k) {
            int j = m_adjncy[k];
            mat[i+ne*j] = 1;
            mat[j+ne*i] = 1;
        }
    }
    fwrite(mat, ne*ne, 1, fmat);
    fclose(fmat);
}


void Mesh3D::write_materials(const std::string& str)
{
    write_materials_v2(str);
}

int Mesh3D::read_materials(const std::string& str)
{
    return read_materials_v2(str);
}

int Mesh3D::read_materials_v1(const std::string& str)
{
    FILE* f = fopen(str.c_str(), "r");
    int nmats;
    char type;
    char *buffer=NULL;
    size_t linesize=0;
    double vs, vp, rho;
    int ngllx, nglly, ngllz;
    double dt, Qp, Qmu;

    getline(&buffer, &linesize, f);
    sscanf(buffer, "%d", &nmats);
    for(int k=0;k<nmats;++k) {
        getline(&buffer, &linesize, f);
        sscanf(buffer, "%c %lf %lf %lf %d %d %d %lf %lf %lf",
               &type, &vp, &vs, &rho, &ngllx, &nglly, &ngllz,
               &dt, &Qp, &Qmu);
        printf("Mat: %2ld : %c vp=%lf vs=%lf Qp=%lf Qmu=%lf\n", m_materials.size(), type, vp, vs, Qp, Qmu);
        m_materials.push_back(Material(type, vp, vs, rho, Qp, Qmu, ngllx));
    }
    free(buffer);
    return nmats;
}

int Mesh3D::read_materials_v2(const std::string& str)
{
    FILE* f = fopen(str.c_str(), "r");
    int nmats;
    char type;
    char *buffer=NULL;
    size_t linesize=0;
    double vs, vp, rho;
    int ngllx;
    double Qp, Qmu;

    getline(&buffer, &linesize, f);
    sscanf(buffer, "%d", &nmats);
    for(int k=0;k<nmats;++k) {
        getline(&buffer, &linesize, f);
        sscanf(buffer, "%c %lf %lf %lf %d %lf %lf",
               &type, &vp, &vs, &rho, &ngllx, &Qp, &Qmu);
        printf("Mat: %2ld : %c vp=%lf vs=%lf\n", m_materials.size(), type, vp, vs);
        m_materials.push_back(Material(type, vp, vs, rho, Qp, Qmu, ngllx));
    }
    free(buffer);
    return nmats;
}

#define TF(e)  (e ? 'T' : 'F')


void Mesh3D::write_materials_v1(const std::string& str)
{
    FILE* f = fopen(str.c_str(), "w");
    int nmats = m_materials.size();
    fprintf(f, "%d\n", nmats);
    for(int k=0;k<nmats;++k) {
        const Material& mat = m_materials[k];
        fprintf(f, "%c %lf %lf %lf %d %d %d %lf %lf %lf\n",
                mat.material_char(),
                mat.Pspeed, mat.Sspeed, mat.rho,
                mat.m_ngll, mat.m_ngll, mat.m_ngll,
                0.0, mat.Qpression, mat.Qmu);
    }
    fprintf(f, "# PML properties\n");
    fprintf(f, "# Filtering? npow,Apow,X?,left?,Y?,Forwrd?,Z?,down?,cutoff freq\n");
    for(int k=0;k<nmats;++k) {
        const Material& mat = m_materials[k];
        if (!mat.is_pml()) continue;
        char px = TF(mat.xwidth!=0.);
        char py = TF(mat.ywidth!=0);
        char pz = TF(mat.zwidth!=0);
        char xl = TF(mat.xwidth<0);
        char yf = TF(mat.ywidth<0);
        char zd = TF(mat.zwidth<0);
        fprintf(f, "F 2 10. %c %c %c %c %c %c 0. 0\n", px, xl, py, yf, pz, zd);
    }
}

void Mesh3D::write_materials_v2(const std::string& str)
{
    FILE* f = fopen(str.c_str(), "w");
    int nmats = m_materials.size();
    fprintf(f, "%d\n", nmats);
    for(int k=0;k<nmats;++k) {
        const Material& mat = m_materials[k];
        fprintf(f, "%c %lf %lf %lf %d %lf %lf\n",
                mat.material_char(),
                mat.Pspeed, mat.Sspeed, mat.rho,
                mat.m_ngll,
                mat.Qpression, mat.Qmu);
    }

    fprintf(f, "# PML properties\n");
    fprintf(f, "# npow,Apow,posX,widthX,posY,widthY,posZ,widthZ,mat\n");
    for(int k=0;k<nmats;++k) {
        const Material& mat = m_materials[k];
        if (!mat.is_pml()) continue;
        fprintf(f, "2 10. %lf %lf %lf %lf %lf %lf %d\n",
                mat.xpos, mat.xwidth,
                mat.ypos, mat.ywidth,
                mat.zpos, mat.zwidth, k);
    }
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
    m_vertex_domains.clear();
    m_vertex_domains.resize(n_vertices(), 0);
    for(int i=0;i<nel;++i) {
        for(int k=m_elems_offs[i];k<m_elems_offs[i+1];++k) {
            int vtx = m_elems[k];
            int domain = m_materials[m_mat[i]].domain();
            m_vertex_to_elem.add_link(vtx, i);
            m_vertex_domains[vtx] |= (1<<domain);
//            printf("VX[%d] dom=%d/%02x, %02x\n", vtx, domain, (int)(1<<domain), m_vertex_domains[vtx]);
        }
    }
//    for(int k=0;k<n_vertices();++k) {
//        printf("VX[%d] dom=%02x\n", k, m_vertex_domains[k]);
//    }
}


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
