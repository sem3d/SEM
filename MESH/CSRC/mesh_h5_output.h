/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#ifndef MESH_H5_OUTPUT_H
#define MESH_H5_OUTPUT_H

class Mesh3D;
class MeshPart;

void output_mesh_part_h5(MeshPart& loc, Mesh3D& mesh);
void output_mesh_part_h5_comm(int part, Mesh3D& mesh);
void output_all_meshes_xmf(int nprocs);
void output_mesh_part_xmf(MeshPart& loc, Mesh3D& mesh);

#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
