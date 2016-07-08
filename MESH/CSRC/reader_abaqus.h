/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// mesh.h : Gestion maillage format SEM
#ifndef _READER_ABAQUS_H_
#define _READER_ABAQUS_H_

#include "mesh.h"
#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include "mesh_common.h"

class MeshReaderAbaqus {
public:
    MeshReaderAbaqus(const char* fname);
    ~MeshReaderAbaqus();
    void parse_file(Mesh3D& mesh);

    bool eof();
protected:
    char curline[2048];
    FILE* m_file;
    int m_line;
    void read_next_line();
    void read_header(Mesh3D& mesh);
    void read_node_section(Mesh3D& mesh);
    void read_elements(Mesh3D& mesh);

    std::map<int,int> m_node_map; // maps file id to mesh ids
    std::vector<std::string> m_elemsets;
};

#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
