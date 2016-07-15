/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include "mesh_h5_output.h"
#include "mesh.h"
#include "h5helper.h"
#include "mesh_common.h"

#include <vector>
#include <cstdio>

using std::vector;


static void output_all_meshes_xmf_part(int nprocs, const char* sfx)
{
    char fname[2048];
    FILE* f;

    snprintf(fname, sizeof(fname), "mesh4spec.%s.xmf", sfx);
    f = fopen(fname,"w");
    fprintf(f, "<?xml version=\"1.0\" ?>\n");
    fprintf(f, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n");
    fprintf(f, "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n");
    fprintf(f, "  <Domain>\n");
    fprintf(f, "    <Grid name=\"mesh\" CollectionType=\"Spatial\" GridType=\"Collection\">\n");
    for(int n=0;n<nprocs;++n)
    {
	fprintf(f, "      <xi:include href=\"mesh4spec.%04d.%s.xmf\" xpointer=\"xpointer(//Xdmf/Domain/Grid)\" />\n", n, sfx);
    }
    fprintf(f, "    </Grid>\n");
    fprintf(f, "  </Domain>\n");
    fprintf(f, "</Xdmf>\n");
    fclose(f);
}

void output_all_meshes_xmf(int nprocs)
{
    output_all_meshes_xmf_part(nprocs, "elems");
    output_all_meshes_xmf_part(nprocs, "faces");
    output_all_meshes_xmf_part(nprocs, "edges");
    output_all_meshes_xmf_part(nprocs, "comms.edges");
    output_all_meshes_xmf_part(nprocs, "comms.faces");
}


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
