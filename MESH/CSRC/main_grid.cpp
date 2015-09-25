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
#include "meshpart.h"
#include "mesh_h5_output.h"
using namespace std;



// Describe directions of presence of PML elements
// When generating a structured axis aligned grid
struct Pml_dirs {
    bool N,S,E,W,U,D;
};



struct RectMesh {
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    double xstep, ystep, zstep;
    Pml_dirs pmls;
};

int pointidx(int i, int j, int k, int nx, int ny, int nz)
{
    return i+j*(nx+1) + k*(nx+1)*(ny+1);
}



void init_rectangular_mesh( RectMesh& desc, Mesh3D& mesh)
{
    assert(desc.xmin<desc.xmax);
    assert(desc.ymin<desc.ymax);
    assert(desc.zmin<desc.zmax);

    if (desc.pmls.U) {} // TODO add some elements along the border

    int nelemx = int( (desc.xmax-desc.xmin)/desc.xstep );
    int nelemy = int( (desc.ymax-desc.ymin)/desc.ystep );
    int nelemz = int( (desc.zmax-desc.zmin)/desc.zstep );

    printf("Creating grid mesh %d x %d x %d\n", nelemx, nelemy, nelemz );
    // Coordinates
    for(int k=0;k<=nelemz;++k) {
	for(int j=0;j<=nelemy;++j) {
	    for(int i=0;i<=nelemx;++i) {
		mesh.add_node( desc.xmin + i*desc.xstep,
			       desc.ymin + j*desc.ystep,
			       desc.zmin + k*desc.zstep );
	    }
	}
    }
    // Elements
    HexElem elem;
    for(int k=0;k<nelemz;++k) {
	for(int j=0;j<nelemy;++j) {
	    for(int i=0;i<nelemx;++i) {
		elem.v[0] = pointidx(i  ,j  ,k  ,nelemx,nelemy,nelemz);
		elem.v[1] = pointidx(i+1,j  ,k  ,nelemx,nelemy,nelemz);
		elem.v[2] = pointidx(i+1,j+1,k  ,nelemx,nelemy,nelemz);
		elem.v[3] = pointidx(i  ,j+1,k  ,nelemx,nelemy,nelemz);
		elem.v[4] = pointidx(i  ,j  ,k+1,nelemx,nelemy,nelemz);
		elem.v[5] = pointidx(i+1,j  ,k+1,nelemx,nelemy,nelemz);
		elem.v[6] = pointidx(i+1,j+1,k+1,nelemx,nelemy,nelemz);
		elem.v[7] = pointidx(i  ,j+1,k+1,nelemx,nelemy,nelemz);
		mesh.add_elem(0, elem);
	    }
	}
    }
}





int main(int argc, char**argv)
{
    Mesh3D mesh;
    RectMesh desc;

    desc.xmin = -100;
    desc.xmax =  100;
    desc.ymin = -100;
    desc.ymax =  100;
    desc.zmin = -100;
    desc.zmax =  100;
    desc.xstep = desc.ystep = desc.zstep = 10;
    desc.pmls.N = false;
    desc.pmls.S = false;
    desc.pmls.E = false;
    desc.pmls.W = false;
    desc.pmls.U = false;
    desc.pmls.D = true;

    int NPROCS=atoi(argv[1]);

    mesh.read_materials("material.input");
    init_rectangular_mesh(desc, mesh);

    mesh.partition_mesh(NPROCS);
    mesh.build_vertex_to_elem_map();

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
