/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// mesh.h : Gestion maillage format SEM
#ifndef _MESH_COMMON_
#define _MESH_COMMON_

#include <cstdio>

// Reads a line skiping it if it starts with #
void getData_line(char **buffer, size_t* linesize, FILE* f)
{
    ssize_t nc;
    nc = getline(buffer, linesize, f);

    for(int k=0;k<100;++k) {
        if(*buffer[0] != '#'){
          break;
        }
        else{
            nc = getline(buffer, linesize, f);
        }
    }
    if (nc<=0 && *buffer) { (*buffer)[0] = 0; }
}


#endif




/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
