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

void output_mesh_part_h5(MeshPart& loc, Mesh3D& mesh)
{
    int i;
    int part = loc.part;
    char fname[2048];
    vector<double> tmpd;
    vector<int> tmpi, tmpi1;

    printf("%04d : number of elements = %d\n", part, loc.n_elems());

    snprintf(fname, sizeof(fname), "mesh4spec.%04d.h5", part);
    hid_t fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    h5h_create_attr(fid, "ndims", 3);
    h5h_create_attr(fid, "n_processors", mesh.n_parts());
    h5h_create_attr(fid, "n_materials", mesh.n_materials());
    h5h_create_attr(fid, "n_elements", loc.n_elems());
    h5h_create_attr(fid, "n_faces", loc.n_faces());
    h5h_create_attr(fid, "n_edges", loc.n_edges());
    h5h_create_attr(fid, "n_vertices", mesh.n_vertices());
    h5h_create_attr(fid, "solid_fluid", false);
    h5h_create_attr(fid, "solid_fluid_loc", false);
    h5h_create_attr(fid, "all_fluid", false);
    h5h_create_attr(fid, "neumann_present", false);
    h5h_create_attr(fid, "neumann_present_loc", false);
    h5h_create_attr(fid, "curve", false);

    mesh.get_local_nodes(loc, tmpd);
    h5h_write_dset(fid, "local_nodes", loc.n_nodes(), 3, &tmpd[0]);

    mesh.get_local_material(loc, tmpi);
    h5h_write_dset(fid, "material", loc.n_elems(), 2, &tmpi[0]);

    mesh.get_local_elements(loc, tmpi);
    h5h_write_dset(fid, "elements", loc.n_elems(), 8, &tmpi[0]);
    //
    h5h_write_dset(fid, "faces", loc.n_elems(), 6, &loc.faces[0]);
    h5h_write_dset(fid, "faces_map", loc.n_elems(), 6, &loc.faces_orient[0]);
    //
    h5h_write_dset(fid, "edges", loc.n_elems(), 12, &loc.edges[0]);
    h5h_write_dset(fid, "edges_map", loc.n_elems(), 12, &loc.edges_orient[0]);
    //
    h5h_write_dset(fid, "vertices", loc.n_elems(), 8, &loc.vertices[0]);
    h5h_write_dset(fid, "vertices_to_global", loc.n_l2g.size(), &loc.n_l2g[0]);
    //
    H5Fclose(fid);
}


void output_mesh_part_h5_comm(int part, Mesh3D& mesh)
{
    int i;
    char fname[2048];
    vector<double> tmpd;
    vector<int> tmpi, tmpi1;

    snprintf(fname, sizeof(fname), "mesh4spec.%04d.h5", part);
    hid_t fid = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

    for(i=0;i<mesh.n_procs;++i) {
	char proc_grp[20];
	snprintf(proc_grp, 20, "Proc%04d", i);
	hid_t grp_id = H5Gcreate(fid, proc_grp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	h5h_create_attr(grp_id, "n_faces", mesh.n_shared_faces());
	h5h_create_attr(grp_id, "n_edges", mesh.n_shared_edges());
	h5h_create_attr(grp_id, "n_vertices", mesh.n_shared_vertices());

	H5Gclose(grp_id);
    }
    H5Fclose(fid);
}

void output_mesh_part_xmf(MeshPart& loc, Mesh3D& mesh)
{
    int i;
    int part = loc.part;
    char fname[2048];
    FILE* f;

    snprintf(fname, sizeof(fname), "mesh4spec.%04d.xmf", part);
    f = fopen(fname,"w");
    fprintf(f, "<?xml version=\"1.0\" ?>\n");
    fprintf(f, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n");
    fprintf(f, "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n");
    fprintf(f, "  <Domain>\n");
    fprintf(f, "    <Grid name=\"mesh.%04d\">\n", part);
    fprintf(f, "      <Topology Type=\"Hexahedron\" NumberOfElements=\"%d\">\n", loc.n_elems());
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%d 8\">\n", loc.n_elems());
    fprintf(f, "mesh4spec.%04d.h5:/elements\n",part);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Topology>\n");
    fprintf(f, "      <Geometry Type=\"XYZ\">\n");
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%d 3\">\n", loc.n_nodes());
    fprintf(f, "mesh4spec.%04d.h5:/local_nodes\n", part);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Geometry>\n");
    fprintf(f, "    </Grid>\n");
    fprintf(f, "  </Domain>\n");
    fprintf(f, "</Xdmf>\n");
    fclose(f);
}

void output_all_meshes_xmf(int nprocs)
{
    int i;
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
