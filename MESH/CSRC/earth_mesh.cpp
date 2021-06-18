/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// earthmesh.cpp : Generates a sphere meshed with hexaedrons and refined layers

#include <cmath>
#include <cassert>
#include "earth_mesh.h"

void EarthMesh::read_params(FILE* inp)
{
    int mat, ndiv;
    double zz;
    char surftype[100];
    char surfname[100];
    char layertype[100];

    fscanf(inp, "%d %lf %d %d", &nlayers, &Z0, &N0, &mat);
    add_layer(2, mat, Z0, N0);
    for(int i=0;i<nlayers;++i) {
        fscanf(inp, "%s", surftype);
        if (surftype[0]=='z') {
            fscanf(inp, "%lf", &zz);
        } else {
            fscanf(inp, "%s", surfname);
        }
        fscanf(inp, "%d %s %d", &mat, layertype, &ndiv);
        // Pour l'instant on ne gere pas raf juste div
        add_layer(0, mat, zz, ndiv);
    }
}

void EarthMesh::add_layer(int type, int mat, double z, int ndiv)
{
    layer_ndiv.push_back(ndiv);
    layer_type.push_back(type);
    layer_z1.push_back(z);
    layer_mat.push_back(mat);
}

void EarthMesh::init_earth()
{
    create_center_cube();

    // Les noeuds du cube centre sont numérotés de 0 à (N0+1)^3-1
    // par I+J*(N0+1)+K*(N0+1)^2
    // On passe N0, DI, DJ à face cube pour qu'il puisse recalculer
    // Les numéros de noeuds de la face I,J,. ou I,.,J ou .,I,J
    // par N0+I*DI+J*DJ
    // Ainsi la premiere face (K=0) est designée par (0,1,N0+1)
    // et indexe (-W,-W,-W) ... (W,W,-W)
    // Mais il faut lui donner une orientation cohérente avec Z<0
    create_face_cube_top();  // Z+
    create_face_cube_bottom(); // Z-
    create_face_cube_left();   // X-
    create_face_cube_right();  // X+
    create_face_cube_front();  // Y-
    create_face_cube_back();   // Y+
}

void EarthMesh::create_center_cube()
{
    Elem el(8);
    points.create_cube(Z0, N0);
    int mat = layer_mat[0];
    for(int k=0;k<=N0;++k) {
        for(int j=0;j<=N0;++j) {
            for(int i=0;i<=N0;++i) {
                int n = i + (N0+1)*(j + (N0+1)*k);
                Point3D& pt = points.pts[n];
                int nn = mesh.add_node(pt.x, pt.y, pt.z);
                assert(n==nn);
            }
        }
    }
    for(int k=0;k<N0;++k) {
        for(int j=0;j<N0;++j) {
            for(int i=0;i<N0;++i) {
                el.v[0] = i+0 + (N0+1)*(j+0 + (N0+1)*(k+0));
                el.v[1] = i+1 + (N0+1)*(j+0 + (N0+1)*(k+0));
                el.v[2] = i+1 + (N0+1)*(j+1 + (N0+1)*(k+0));
                el.v[3] = i+0 + (N0+1)*(j+1 + (N0+1)*(k+0));
                el.v[4] = i+0 + (N0+1)*(j+0 + (N0+1)*(k+1));
                el.v[5] = i+1 + (N0+1)*(j+0 + (N0+1)*(k+1));
                el.v[6] = i+1 + (N0+1)*(j+1 + (N0+1)*(k+1));
                el.v[7] = i+0 + (N0+1)*(j+1 + (N0+1)*(k+1));
                mesh.add_elem(mat, el);
            }
        }
    }
}

void EarthMesh::create_face_cube_top()
{
}
void EarthMesh::create_face_cube_bottom()
{
}
void EarthMesh::create_face_cube_left()
{
}
void EarthMesh::create_face_cube_right()
{
}
void EarthMesh::create_face_cube_front()
{
}
void EarthMesh::create_face_cube_back()
{
}


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
