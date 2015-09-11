/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include "mesh_h5_output.h"
#include "mesh.h"
#include "h5helper.h"

#include <vector>
#include <cstdio>

using std::vector;


void output_all_meshes_xmf(int nprocs)
{
    char fname[2048];
    FILE* f;

    snprintf(fname, sizeof(fname), "mesh4spec.xmf");
    f = fopen(fname,"w");
    fprintf(f, "<?xml version=\"1.0\" ?>\n");
    fprintf(f, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n");
    fprintf(f, "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n");
    fprintf(f, "  <Domain>\n");
    fprintf(f, "    <Grid name=\"mesh\" CollectionType=\"Spatial\" GridType=\"Collection\">\n");
    for(int n=0;n<nprocs;++n)
    {
	fprintf(f, "      <xi:include href=\"mesh4spec.%04d.xmf\" xpointer=\"xpointer(//Xdmf/Domain/Grid)\" />\n", n);
    }
    fprintf(f, "    </Grid>\n");
    fprintf(f, "  </Domain>\n");
    fprintf(f, "</Xdmf>\n");
    fclose(f);
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
