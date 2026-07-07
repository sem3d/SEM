/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

//mat.dat
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include "material.h"
#include "mesh.h"
#include "meshpart.h"
#include "mesh_h5_output.h"
#include "mesh_grid.h"
#include "reader_abaqus.h"
#include "reader_ideas.h"
#include "mesh_common.h"
#include "read_input.h"
#include "sem_input.h"
#include "earth_mesh.h"
#include "mesh_pml_extrude.h"

static bool mesh_has_pml_materials(const Mesh3D& mesh)
{
    for (size_t k=0;k<mesh.m_materials.size();++k)
        if (is_dm_pml(mesh.m_materials[k].domain())) return true;
    return false;
}

// For an imported mesh whose PML materials are already present (declared P/L in mater.in
// but without descriptors), fill each PML material's pos/width and associated material from
// the geometry: pos/width come from how far the PML material's elements extend beyond the
// interior (non-PML) domain on each axis; the associated interior material is inferred from
// a shared node.
static void derive_pml_descriptors_from_geometry(Mesh3D& mesh)
{
    size_t nmat = mesh.m_materials.size();
    size_t ne = mesh.n_elems();
    const double BIG = 1e300;
    std::vector<double> mnx(nmat,BIG), mxx(nmat,-BIG), mny(nmat,BIG), mxy(nmat,-BIG), mnz(nmat,BIG), mxz(nmat,-BIG);
    double ixmn=BIG,ixmx=-BIG,iymn=BIG,iymx=-BIG,izmn=BIG,izmx=-BIG;
    std::map<index_t,int> node_interior_mat; // node -> a non-PML material (for assoc)

    for (size_t e=0;e<ne;++e) {
        int mat = mesh.m_mat[e];
        bool pml = is_dm_pml(mesh.m_materials[mat].domain());
        index_t nodes[8]; mesh.get_elem_nodes((index_t)e, nodes);
        for (int i=0;i<8;++i) {
            double x=mesh.m_xco[nodes[i]], y=mesh.m_yco[nodes[i]], z=mesh.m_zco[nodes[i]];
            if (x<mnx[mat])mnx[mat]=x; if (x>mxx[mat])mxx[mat]=x;
            if (y<mny[mat])mny[mat]=y; if (y>mxy[mat])mxy[mat]=y;
            if (z<mnz[mat])mnz[mat]=z; if (z>mxz[mat])mxz[mat]=z;
            if (!pml) {
                if (x<ixmn)ixmn=x; if (x>ixmx)ixmx=x;
                if (y<iymn)iymn=y; if (y>iymx)iymx=y;
                if (z<izmn)izmn=z; if (z>izmx)izmx=z;
                node_interior_mat[nodes[i]] = mat;
            }
        }
    }
    if (ixmx<ixmn) { printf("ERR: mesh has PML materials but no interior (non-PML) elements\n"); exit(1); }
    double tx=1e-6*(ixmx-ixmn), ty=1e-6*(iymx-iymn), tz=1e-6*(izmx-izmn);

    for (size_t k=0;k<nmat;++k) {
        if (!is_dm_pml(mesh.m_materials[k].domain())) continue;
        if (mxx[k]<mnx[k]) continue; // material unused in the mesh
        double xp=0,xw=0,yp=0,yw=0,zp=0,zw=0;
        if      (mnx[k] < ixmn-tx) { xp=ixmn; xw=mnx[k]-ixmn; }
        else if (mxx[k] > ixmx+tx) { xp=ixmx; xw=mxx[k]-ixmx; }
        if      (mny[k] < iymn-ty) { yp=iymn; yw=mny[k]-iymn; }
        else if (mxy[k] > iymx+ty) { yp=iymx; yw=mxy[k]-iymx; }
        if      (mnz[k] < izmn-tz) { zp=izmn; zw=mnz[k]-izmn; }
        else if (mxz[k] > izmx+tz) { zp=izmx; zw=mxz[k]-izmx; }
        mesh.m_materials[k].set_pml_borders(xp,xw,yp,yw,zp,zw);

        int assoc = -1;
        for (size_t e=0;e<ne && assoc<0;++e) {
            if (mesh.m_mat[e]!=(int)k) continue;
            index_t nodes[8]; mesh.get_elem_nodes((index_t)e, nodes);
            for (int i=0;i<8;++i) {
                std::map<index_t,int>::iterator it=node_interior_mat.find(nodes[i]);
                if (it!=node_interior_mat.end()) { assoc=it->second; break; }
            }
        }
        mesh.m_materials[k].associated_material = assoc;
        printf("PML material %zu (from geometry): x(pos=%g,w=%g) y(pos=%g,w=%g) z(pos=%g,w=%g) assoc=%d\n",
               k, xp,xw,yp,yw,zp,zw, assoc);
        if (assoc<0) printf("  WARNING: could not infer associated interior material for PML %zu\n", k);
    }
}

void handle_on_the_fly(Mesh3D& mesh)
{
    FILE* fparams;
    RectMesh desc;
    fparams = fopen("mat.dat","r");
    desc.read_params_old(fparams);
    fclose(fparams);
    desc.init_rectangular_mesh(mesh);
}



void handle_ideas_file(Mesh3D& mesh)
{
    int numfiles;
    char fname[2048];

    printf("\nHow many files ?\n");
    scanf("%d", &numfiles);
    for(int k=0;k<numfiles;++k) {
        printf("File %d name ?\n", k+1);
        scanf("%2000s", fname);
        MeshReaderIdeas  reader(fname);
        reader.parse_file(mesh, fname);
    }
}

void handle_abaqus_file(Mesh3D& mesh)
{
    int numfiles;
    char fname[2048];

    printf("\nHow many files ?\n");
    scanf("%d", &numfiles);
    for(int k=0;k<numfiles;++k) {
        printf("File %d name ?\n", k+1);
        scanf("%2000s", fname);
        MeshReaderAbaqus  reader(fname);
        reader.parse_file(mesh);
    }
}

void handle_hdf5_file(Mesh3D& mesh)
{
    int numfiles;
    char fname[2048];

     printf("\nHow many files ?\n");
     scanf("%d", &numfiles);
     for(int k=0;k<numfiles;++k) {
         printf("File %d name ?\n", k+1);
         scanf("%2000s", fname);
	 mesh.read_mesh_file(fname);
	 }
}

void handle_earth_chunk(Mesh3D& mesh)
{
}

void handle_full_earth(Mesh3D& mesh)
{
    /* earth.dat format:
NLAYERS  Z0 N0 MAT0
z Z1 mat div NDIV
z Z2 mat raf 1
z Z3 mat div 5
s surfname mat div 6

z: indicate a fixed scalar value
s: interpolate surface depth from surfname=h5file/dataset 
mat: material number
N0 : number of cell of center cube across one edge
MAT0: material of the center cube
div/raf either divide the layer or refine
Z is distance from center and should appear in increasing order
Z0<Z1<Z2<Z3<surf

The center cube has Nc= Ninit*Ninit*Ninit cells.
    */
    FILE* fparams;
    EarthMesh earth(mesh);
    fparams = fopen("earth.dat","r");
    earth.read_params(fparams);
    fclose(fparams);
    earth.init_earth();
}

/// Emulates old mesher interface

int main(int argc, char**argv)
{
    Mesh3D mesh;
    int NPROCS;
    int choice;
    char *buffer=NULL;
    //FILE* f = fopen("mesh.input", "r");
    //char *buffer=NULL;
    size_t linesize=0;


    mesh.debug = true;
    printf("-------------------------------------------------\n");
    printf("-------------------------------------------------\n");
    printf("-----                                       -----\n");
    printf("----- Construction of input files for SEM3D -----\n");
    printf("-----                                       -----\n");
    printf("-------------------------------------------------\n");
    printf("-------------------------------------------------\n");
    if (mesh.debug) {
        printf("\n    DEBUG MODE    \n\n");
    }
    printf("\n   --> How many procs for the run ?\n");

    getData_line(&buffer, &linesize, stdin);

    sscanf(buffer,"%d", &NPROCS);

    printf("             %d processor(s)\n", NPROCS);
    printf(" \n\n");
    printf("  --> Which Initial Mesh?\n");
    printf("      1- On the fly\n");
    printf("      2- Abaqus from Cubit\n");
    printf("      3- Ideas (.unv) files\n");
    printf("      4- HDF5 Hex8 files\n");
    printf("      5- Earth Chunk\n");
    printf("      6- Full earth\n");

    getData_line(&buffer, &linesize, stdin);

    sscanf(buffer,"%d", &choice);
    printf("            Your choice is %d \n", choice);
    printf(" \n\n");

    sem_config_t config;
    int err;
    read_sem_config(&config, 0, 3, "input.spec", &err);
    dump_config(&config);

    switch(choice) {
    case 1:
        mesh.read_materials("mater.in");
        handle_on_the_fly(mesh);
        mesh.write_materials("material.input");
        break;
    case 2:
        if (access("mater.in", F_OK)==0) mesh.read_materials("mater.in", false);
        else mesh.read_materials("material.input");
        handle_abaqus_file(mesh);
        break;
    case 3:
        if (access("mater.in", F_OK)==0) mesh.read_materials("mater.in", false);
        else mesh.read_materials("material.input");
        handle_ideas_file(mesh);
        break;
    case 4:
        if (access("mater.in", F_OK)==0) mesh.read_materials("mater.in", false);
        else mesh.read_materials("material.input");
        handle_hdf5_file(mesh);
        break;
    case 5:
        handle_earth_chunk(mesh);
        break;
    case 6:
        handle_full_earth(mesh);
        break;
    default:
        break;
    };

    // Anisotropy (Fluid_Aniso/Cstar_Fluid deftype, like Hooke_Aniso for solids) is a
    // material.spec-only, runtime concern: fluid materials stay tagged 'F'/DM_FLUID_CG in
    // the mesh and material.input, exactly like anisotropic solids stay 'S'/DM_SOLID_CG.
    // No mesh-side domain upgrade needed (there used to be one baking DM_FLUID_CG_ANISO
    // here; that domain no longer exists in the solver -- aniso is dom%aniso on the
    // regular fluid domain now).

    // PML source selection for imported meshes (cases 2/3/4):
    //  - if PML materials are already defined (mater.in flagged them P/L, so the imported
    //    mesh already contains the PML elements) -> derive their descriptors from geometry
    //    and ignore pml.input;
    //  - otherwise, if pml.input exists -> extrude PML layers.
    // Only material.input is (re)written; PMLs are standard isotropic P/L materials there,
    // so material.spec (if any) only needs to define the interior materials.
    // Case 1 (on the fly) keeps its PML from mat.dat; pml.input is ignored there.
    if (choice==2 || choice==3 || choice==4) {
        if (mesh_has_pml_materials(mesh)) {
            derive_pml_descriptors_from_geometry(mesh);
            if (access("pml.input", F_OK)==0)
                printf("WARNING: pml.input ignored (the mesh already declares PML materials in mater.in)\n");
        } else {
            PmlSpec pmlspec;
            if (read_pml_input("pml.input", pmlspec) && pmlspec.any())
                extrude_pml(mesh, pmlspec);
        }
        mesh.write_materials("material.input");
    } else if (choice==1 && access("pml.input", F_OK)==0) {
        printf("WARNING: pml.input ignored for 'on the fly' meshes (PML comes from mat.dat)\n");
    }

    //mesh.write_materials("material.input");
    mesh.define_associated_materials();

    mesh.generate_output(NPROCS, &config);
    return 0;
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
