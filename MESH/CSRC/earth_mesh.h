/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// earthmesh.h : Generates a sphere meshed with hexaedrons and refined layers
#ifndef _EARTHMESH_H_
#define _EARTHMESH_H_

#include <cstdio>
#include <vector>
#include "mesh.h"
#include "point3d.h"

class EarthMesh {
public:
    void hop();
    EarthMesh(Mesh3D& mesh_):
        mesh(mesh_)
        {}

    void read_params(FILE* inp);
    void init_earth();
    void add_layer(int type, int mat, double z, int ndiv);

    void create_center_cube();
    void create_face_cube_top();  // Z+
    void create_face_cube_bottom(); // Z-
    void create_face_cube_left();   // X-
    void create_face_cube_right();  // X+
    void create_face_cube_front();  // Y-
    void create_face_cube_back();   // Y+

    int N0;
    int nlayers;
    double Z0;
    // Those vectors all have nlayers entries
    // First layer is cube
    std::vector<int>    layer_ndiv;  // Number of division for that layer (if type div, else 0)
    std::vector<int>    layer_type;  // Type of layer div(0) raf(1) centercube(2)
    std::vector<double> layer_z1;    // Top of layer
    std::vector<int>    layer_mat;   // Material number of layer

    Mesh3D& mesh;
    PointsArray_3D  points;
};


#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
