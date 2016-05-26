/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include <cstdio>
#include <cstdlib>
#include "material.h"
#include "mesh.h"
#include "meshpart.h"
#include "mesh_h5_output.h"
#include "mesh_grid.h"
#include "reader_abaqus.h"
#include "reader_ideas.h"

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
}

void handle_earth_chunk(Mesh3D& mesh)
{
}

/// Emulates old mesher interface

int main(int argc, char**argv)
{
    Mesh3D mesh;
    int NPROCS;
    int choice;

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
    scanf("%d", &NPROCS);
    printf("             Your choice is %d \n", NPROCS);

    printf(" ****************************************\n");
    printf(" ****************************************\n");
    printf("  WHICH INITIAL MESH?\n");
    printf("      1- On the fly\n");
    printf("      2- Abaqus from Cubit\n");
    printf("      3- Ideas (.unv) files\n");
    printf("      4- HDF5 Hex8 files\n");
    printf("      5- Earth Chunk\n");
    scanf("%d", &choice);
    printf("            Your choice is %d \n", choice);


    switch(choice) {
    case 1:
        mesh.read_materials("mater.in");
        handle_on_the_fly(mesh);
        break;
    case 2:
        mesh.read_materials("mater.in");
        handle_abaqus_file(mesh);
        break;
    case 3:
    	mesh.read_materials("mater.in");
        handle_ideas_file(mesh);
        break;
    case 4:
        handle_hdf5_file(mesh);
        break;
    case 5:
        handle_earth_chunk(mesh);
        break;
    default:
        break;
    };

    mesh.write_materials("material.input");

    mesh.generate_output(NPROCS);
    return 0;
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
