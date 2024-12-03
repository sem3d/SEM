/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#ifndef MESH_GRID_H
#define MESH_GRID_H

#include <cstdio>
#include <cstdint>
#include "meshbase.h"

// Describe directions of presence of PML elements
// When generating a structured axis aligned grid
struct Pml_dirs {
    bool N,S,E,W,U,D;
};

struct mirrors_dirs {
    bool W, E, S, N, D, U;
    int nW, nE, nS, nN, nD, nU;
};

class Mesh3D;

struct RectMesh {
    double xmin, xmax;
    double ymin, ymax;
    double zmax, zmin;
    double xmin0, xmax0;
    double ymin0, ymax0;
    double zmax0, zmin0;
    double xstep, ystep;
    int has_pml;
    Pml_dirs pmls;

    int has_mirrors;
    mirrors_dirs mrrs;

    int has_sph;

    int nlayers;
    int npml; // Number of layers of pmls
    double zheight;
    double *thickness;
    int    *nsteps;
    //int ngll_pml;
    int elem_shape;

    void read_params(FILE* fparam);
    void read_params_old(FILE* fparam);
    void init_rectangular_mesh(Mesh3D& mesh);
    void apply_pml_borders();
protected:

    int nelemx, nelemy, nelemz;
    int pointidx(int i, int j, int k);
    int pointidx27(int i, int j, int k);
    int get_mat(Mesh3D& mesh, int layer, bool W, bool E, bool S, bool N, bool U, bool D);
    void emit_free_face(Surface* surf, int dom, const Elem& elem,
                        bool W, bool E, bool S, bool N, bool U, bool D);
    void emit_free_face(Surface* surf, int dom, const Elem& elem, int facenum);
    void emit_mirr_face(Surface* surf, int dom, const Elem& elem,
                        bool W, bool E, bool S, bool N, bool U, bool D);
    void emit_mirr_face(Surface* surf, int dom, const Elem& elem, int facenum, bool flp);
    int create_linear_grid_nodes(Mesh3D& mesh);
    int create_quadratic_grid_nodes(Mesh3D& mesh);
    int create_linear_grid_nodes_sph(Mesh3D& mesh); // spherical shell transformation

    void create_linear_element(Elem& elem, int i, int j, int k);
    void create_quadratic_element(Elem& elem, int i, int j, int k);
};

#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
