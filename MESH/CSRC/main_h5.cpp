/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include "material.h"
#include "mesh.h"
#include "mesh_h5_output.h"
#include "meshpart.h"

using namespace std;



int main(int argc, char**argv)
{
    Mesh3D mesh;

    if (argc!=3) {
      printf("Usage:\n");
      printf("sem_part_h5 NPROCS  meshfile.h5\n");
      printf("input.spec and material.input file need to be in current directory\n");
      exit(1);
    }
    int NPROCS=atoi(argv[1]);

    mesh.read_materials("material.input");
    mesh.read_mesh_file(argv[2]);
    mesh.build_vertex_to_elem_map();
    mesh.partition_mesh(NPROCS);

    for(int part=0;part<NPROCS;++part) {
	Mesh3DPart loc(mesh, part);

	loc.compute_part();
	loc.output_mesh_part();
	loc.output_mesh_part_xmf();
    }
    output_all_meshes_xmf(NPROCS);
    return 0;
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
