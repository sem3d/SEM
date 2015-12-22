/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// mesh.h : Gestion maillage format SEM
#ifndef _READER_IDEAS_H_
#define _READER_IDEAS_H_

#include "mesh.h"
#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include "../../COMMON/read_unv.hpp"


class MeshReaderIdeas {
public:
    MeshReaderIdeas(const char* fname);
    ~MeshReaderIdeas();
    void parse_file(Mesh3D& mesh, const char* fname);

    bool eof();
protected:
    char curline[2048];
    FILE* m_file;
    int m_line;
    lsnodes  nodes;
	lselems  elems;

    void read_node_section(Mesh3D& mesh, lsnodes& nodes);
    void read_elements(Mesh3D& mesh, lselems& elems);

    std::map<int,int> m_node_map; // maps file id to mesh ids
    std::vector<int> m_elemsetsID;

};

#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
