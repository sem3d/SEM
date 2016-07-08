/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */


#include <cstdio>
#include <cstring>

#include <iostream>      // cout, cerr
#include <fstream>       // ifstream
#include <sstream>       // stringstream
#include <unordered_map> // unordered_map
#include <algorithm>     // find

#include "mesh.h"
#include "reader_ideas.h"
#include "mesh_common.h"

using namespace std;

MeshReaderIdeas::MeshReaderIdeas(const char* fname)
{
    m_file =  fopen(fname,"r");
    m_line = 0;
}

MeshReaderIdeas::~MeshReaderIdeas()
{
    fclose(m_file);
}

void MeshReaderIdeas::parse_file(Mesh3D& mesh, const char* fname)
{
	read_unv_mesh(fname, nodes, elems);
	read_node_section(mesh, nodes);
	read_elements(mesh, elems);
}

bool MeshReaderIdeas::eof()
{
    return feof(m_file)!=0;
}

void MeshReaderIdeas::read_node_section(Mesh3D& mesh, lsnodes& nodes)
{
    long nn=0;
    int  nid;

	for(unsigned int i = 0; i < nodes.size(); ++i)
	{
		nid = mesh.add_node(get<0>(nodes[i]), get<1>(nodes[i]), get<2>(nodes[i]));
		m_node_map[i] = nid;
		++nn;
	}
}

void MeshReaderIdeas::read_elements(Mesh3D& mesh, lselems& elems)
{
    int info[5];
    int nNodes;
    int num;
    int n[8];
    HexElem el;
    int mat_id;

    long nn=0;

	for(unsigned int i = 0; i < elems.size(); ++i)
	{
		for(int k=0;k<8;++k) {
		   el.v[k] = get<2>(elems[i])[k]	;
		}
		mat_id = get<1>(get<3>(elems[i]));
		mesh.add_elem(mat_id, el);
		++nn;
	}
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
