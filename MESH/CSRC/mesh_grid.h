/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#ifndef MESH_GRID_H
#define MESH_GRID_H

#include <cstdio>

// Describe directions of presence of PML elements
// When generating a structured axis aligned grid
struct Pml_dirs {
    bool N,S,E,W,U,D;
};

class Mesh3D;

struct RectMesh {
    double xmin, xmax;
    double ymin, ymax;
    double zmax;
    double xstep, ystep;
    int has_pml;
    Pml_dirs pmls;

    int nlayers;
    double zheight;
    double *thickness;
    int    *nsteps;
    int ngll_pml;
    int elem_shape;

    void read_params(FILE* fparam);
    void read_params_old(FILE* fparam);
    void init_rectangular_mesh(Mesh3D& mesh);
    void apply_pml_borders();
protected:

    int nelemx, nelemy, nelemz;
    int pointidx(int i, int j, int k);
    int get_mat(Mesh3D& mesh, int layer, bool W, bool E, bool S, bool N, bool U, bool D);
};

#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
