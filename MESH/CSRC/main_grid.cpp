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




int main(int argc, char**argv)
{
    Mesh3D mesh;
    RectMesh desc;
    FILE* fparams;

    if (argc<2) {
        printf("Usage: part_sem_grid NPROCS [param.dat]\n");
        exit(1);
    }
    int NPROCS=atoi(argv[1]);
    if (argc>2) {
        fparams = fopen(argv[2],"r");
        desc.read_params_old(fparams);
        fclose(fparams);
    } else {
        desc.read_params_old(stdin);
    }

    mesh.read_materials("mater.in");
    desc.init_rectangular_mesh(mesh);
    mesh.write_materials("material.input");
    mesh.build_vertex_to_elem_map();
    mesh.partition_mesh(NPROCS);

    for(int part=0;part<NPROCS;++part) {
	Mesh3DPart loc(mesh, part);

	loc.compute_part();
	loc.output_mesh_part();
	loc.output_mesh_part_xmf();
    }
    output_all_meshes_xmf(NPROCS);
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
