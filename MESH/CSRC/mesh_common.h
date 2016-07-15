/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// mesh.h : Gestion maillage format SEM
#ifndef _MESH_COMMON_
#define _MESH_COMMON_

#include <vector>
#include <string>
#include <cstdio>
/*
#include <map>
#include "material.h"
#include "h5helper.h"
#include "vertex_elem_map.h"
#include "meshbase.h"
#include "aabb.h"
*/

namespace mesh_common{
	size_t getData_line(char **buffer, size_t linesize, FILE* f);
}

#endif

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
