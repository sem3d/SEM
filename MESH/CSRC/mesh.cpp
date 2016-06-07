/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include <cstdio>
#include "mesh.h"
#include "metis.h"
#include <map>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include "mesh_h5_output.h"
#include "meshpart.h"

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

int Mesh3D::add_elem(int mat_idx, const Elem& el)
{
    // Builds elem<->vertex graph
    for(int i=0;i<el.N;++i) {
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
    int corrMod;
	double corrL_x;
	double corrL_y;
	double corrL_z;
	int rho_margiF;
	double rho_var;
	int lambda_margiF;
	double lambda_var;
	int mu_margiF;
	double mu_var;
	int seedStart;

    getline(&buffer, &linesize, f);
    sscanf(buffer, "%d", &nmats);
    for(int k=0;k<nmats;++k) {
        getline(&buffer, &linesize, f);
        sscanf(buffer, "%c %lf %lf %lf %d %lf %lf",
               &type, &vp, &vs, &rho, &ngllx, &Qp, &Qmu);
        printf("Mat: %2ld : %c vp=%lf vs=%lf\n", m_materials.size(), type, vp, vs);
        if(strcmp(&type,"R")){
        	getline(&buffer, &linesize, f);
        	sscanf(buffer, "%d %lf %lf %lf %d %lf %d %lf %d %lf %d",
        	               &corrMod,
						   &corrL_x, &corrL_y, &corrL_z,
						   &rho_margiF, &rho_var,
						   &lambda_margiF, &lambda_var,
						   &mu_margiF, &mu_var,
						   &seedStart);
        	printf("     corrMod = %d, cL_x = %lf, cL_y = %lf, cL_z = %lf, seedStart = %d\n",
        			corrMod, corrL_x, corrL_y, corrL_z, seedStart);
        	m_materials.push_back(Material(type, vp, vs, rho, Qp, Qmu, ngllx,
        			    corrMod,
        				corrL_x, corrL_y, corrL_z,
        				rho_margiF, rho_var,
        				lambda_margiF, lambda_var,
        				mu_margiF, mu_var,
						seedStart));
        }
        else{
        	m_materials.push_back(Material(type, vp, vs, rho, Qp, Qmu, ngllx));
        }
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
        		mat.cinitial_type,
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
                mat.zpos, mat.zwidth, mat.associated_material);
    }
    fprintf(f, "# Random properties\n");
    fprintf(f, "# corrMod, corrL_x, corrL_y, corrL_z, rho_margiF, rho_CV, kappa_margiF, kappa_CV, mu_margiF, mu_CV, seedStart\n");
    for(int k=0;k<nmats;++k) {
        const Material& mat = m_materials[k];
        if (strcmp(&mat.cinitial_type,"R") != 0) continue;
        fprintf(f, "%d %lf %lf %lf %d %lf %d %lf %d %lf %d\n",
                mat.corrMod,
				mat.corrL_x, mat.corrL_y, mat.corrL_z,
                mat.rho_margiF, mat.rho_var,
				mat.lambda_margiF, mat.lambda_var,
				mat.mu_margiF, mat.mu_var, mat.seedStart);
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
    set_control_nodes(8);
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
    set_control_nodes(27);
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
            int mat = m_mat[i];
            int domain = m_materials[mat].domain();
            m_vertex_to_elem.add_link(vtx, i);
            m_vertex_domains[vtx] |= (1<<domain);

            // Update bounding boxen
            m_bbox[mat].update_bounds(Vec3(m_xco[vtx],m_yco[vtx],m_zco[vtx]));
//            printf("VX[%d] dom=%d/%02x, %02x\n", vtx, domain, (int)(1<<domain), m_vertex_domains[vtx]);
        }
    }
//    for(int k=0;k<n_vertices();++k) {
//        printf("VX[%d] dom=%02x\n", k, m_vertex_domains[k]);
//    }
}

void Mesh3D::build_sf_interface()
{
    PFace fc;
    Surface* sf = get_surface("sf");
    for(int el=0;el<n_elems();++el) {
        int dom0 = get_elem_domain(el);
        for(int k=m_xadj[el];k<m_xadj[el+1];++k) {
            int neighbour = m_adjncy[k];
            int dom1 = get_elem_domain(neighbour);
            // Make sure the face normal points inside fluid or fluidpml domain
            if ((dom0==DM_FLUID && dom1==DM_SOLID) ||
                (dom0==DM_FLUID_PML && dom1==DM_SOLID_PML))
            {
                if (get_common_face(el, neighbour, fc)) {
                    fc.set_domain(dom0);
                    sf->add_face(fc,0);
                    fc.set_domain(dom1);
                    fc.orient = -fc.orient;
                    sf->add_face(fc,0);
                }
            }
            if ((dom1==DM_FLUID && dom0==DM_SOLID) ||
                (dom1==DM_FLUID_PML && dom0==DM_SOLID_PML))
            {
                if (get_common_face(neighbour, el, fc)) {
                    fc.set_domain(dom1);
                    sf->add_face(fc,0);
                    fc.set_domain(dom0);
                    fc.orient = -fc.orient;
                    sf->add_face(fc,0);
                }
            }
        }
    }
}

bool Mesh3D::get_common_face(int e0, int e1, PFace& fc)
{
    int nodes0[8];
    int nodes1[8];
    int face0[4];
    std::set<int> inter;
    get_elem_nodes(e0, nodes0);
    get_elem_nodes(e1, nodes1);
    std::set<int> snodes1(nodes1,nodes1+8);
    // walks faces from e0 and return fc if found in e1 with orientation from e0
    for(int nf=0;nf<6;++nf) {
        for(int p=0;p<4;++p) {
            face0[p] = nodes0[RefFace[nf].v[p]];
        }
        std::set<int> sface0(face0,face0+4);
        inter.clear();
        std::set_intersection(sface0.begin(),sface0.end(),snodes1.begin(),snodes1.end(),
                              std::inserter(inter,inter.begin()));
        if (inter.size()==4) {
            fc.set_face(face0);
            return true;
        }
    }
    return false;
}

void Mesh3D::get_neighbour_elements(int nn, const int* n, std::set<int>& elemset) const
{
    std::set<int> elems, temp;
    m_vertex_to_elem.vertex_to_elements(n[0], elemset);
    for(int k=1;k<nn;++k) {
        int vertex_id = n[k];
        elems.clear();
        temp.clear();
        m_vertex_to_elem.vertex_to_elements(n[k], elems);
        std::set_intersection(elemset.begin(), elemset.end(),
                              elems.begin(), elems.end(),
                              std::inserter(temp, temp.begin()));
        elemset.swap(temp);
    }
}

void Mesh3D::save_bbox()
{
    FILE* fbbox;
    fbbox = fopen("domains.txt", "w");
    map<int,AABB>::const_iterator bbox;
    for(bbox=m_bbox.begin();bbox!=m_bbox.end();++bbox) {
        //fprintf(fbbox, "%3d %8.3g %8.3g %8.3g %8.3g %8.3g %8.3g\n", bbox->first,
        fprintf(fbbox, "%3d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", bbox->first,
                bbox->second.min[0],
                bbox->second.min[1],
                bbox->second.min[2],
                bbox->second.max[0],
                bbox->second.max[1],
                bbox->second.max[2]);
    }
}

void Mesh3D::generate_output(int nprocs)
{
    build_vertex_to_elem_map();
    save_bbox();
    partition_mesh(nprocs);
    // partition_mesh builds adjacency map that is used by build_sf_interface
    // Later on, we will want to treat SF interfaces like all others.
    // That is when we will be able to generate normals for all surfaces
    // build_sf_interface();

    for(int part=0;part<nprocs;++part) {
	Mesh3DPart loc(*this, part);

	loc.compute_part();
	loc.output_mesh_part();
	loc.output_mesh_part_xmf();
    }
    output_all_meshes_xmf(nprocs);
}


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
