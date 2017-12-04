/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include <cstdio>
#include "mesh.h"
#include "reader_abaqus.h"
#include <cstring>

MeshReaderAbaqus::MeshReaderAbaqus(const char* fname)
{
    m_file =  fopen(fname,"r");
    m_line = 0;
}

MeshReaderAbaqus::~MeshReaderAbaqus()
{
    fclose(m_file);
}

void MeshReaderAbaqus::parse_file(Mesh3D& mesh)
{
    read_next_line();
    do {
        if (strncmp(curline, "*HEADING", 8)==0) {
            read_header(mesh);
        } else if (strncmp(curline, "*NODE", 5)==0) {
            read_node_section(mesh);
        } else if (strncmp(curline, "*ELEMENT", 8)==0) {
            read_elements(mesh);
        } else if (strncmp(curline, "*PART", 5)==0) {
            read_next_line(); //read_part(mesh);
        } else {
            read_next_line();
        }
    } while (!eof());
}

bool MeshReaderAbaqus::eof()
{
    return feof(m_file)!=0;
}

void MeshReaderAbaqus::read_header(Mesh3D& mesh)
{
    do {
        read_next_line();
        if (curline[0]=='*') break;
    } while(!eof());
}

void MeshReaderAbaqus::read_next_line()
{
    fgets(curline, 2048, m_file);
    m_line++;
    //printf(curline);
}

void MeshReaderAbaqus::read_node_section(Mesh3D& mesh)
{
    int num;
    double x, y, z;
    long nn=0;
    do {
        read_next_line();
        if (curline[0]=='*') break;
        sscanf(curline, "%d, %lf, %lf, %lf", &num, &x, &y, &z);
        int nid = mesh.add_node(x, y, z);
        m_node_map[num] = nid;
        ++nn;
    } while(1);
    printf("Read %ld nodes (line=%d [%s])\n", nn, m_line, curline);
}

void MeshReaderAbaqus::read_elements(Mesh3D& mesh)
{
    char type[30];
    char elset[256];

    sscanf(curline, "*ELEMENT, TYPE=%[^,], ELSET=%s", &type, &elset);
    if (strcmp(type, "C3D8R")!=0) {
        printf("Element type %s ignored\n", type);
        do {
            read_next_line();
            if (eof()) break;
        } while(curline[0]!='*');
        return;
    }
    mesh.set_control_nodes(8);
    int mat_id = m_elemsets.size();
    m_elemsets.push_back(elset);
    long nn=0;
    do {
        int num;
        int n[8];
        HexElem el;
        read_next_line();
        if (curline[0]=='*') break;
        if (eof()) break;
        sscanf(curline, "%d, %d, %d, %d, %d, %d, %d, %d, %d", &num,
              &n[0], &n[1], &n[2], &n[3], &n[4], &n[5], &n[6], &n[7]);
        for(int k=0;k<8;++k) {
            el.v[k] = m_node_map[n[k]];
        }
        mesh.add_elem(mat_id, el);
        ++nn;
    } while(1);
    printf("Read %ld elements from %s (line=%d)\n", nn, elset, m_line);
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
