
#include <vector>
#include <cassert>
#include <cstdio>
#include "material.h"
#include "mesh.h"
#include "h5helper.h"
using namespace std;



// Describe directions of presence of PML elements
// When generating a structured axis aligned grid
struct Pml_dirs {
    bool N,S,E,W,U,D;
};



struct RectMesh {
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    double xstep, ystep, zstep;
    Pml_dirs pmls;
};

int pointidx(int i, int j, int k, int nx, int ny, int nz)
{
    return i+j*(nx+1) + k*(nx+1)*(ny+1);
}



void init_rectangular_mesh( RectMesh& desc, Mesh3D& mesh)
{
    assert(desc.xmin<desc.xmax);
    assert(desc.ymin<desc.ymax);
    assert(desc.zmin<desc.zmax);

    if (desc.pmls.U) {} // TODO add some elements along the border

    int nelemx = int( (desc.xmax-desc.xmin)/desc.xstep );
    int nelemy = int( (desc.ymax-desc.ymin)/desc.ystep );
    int nelemz = int( (desc.zmax-desc.zmin)/desc.zstep );

    // Coordinates
    for(int k=0;k<=nelemz;++k) {
	for(int j=0;j<=nelemz;++j) {
	    for(int i=0;i<=nelemz;++i) {
		mesh.add_node( desc.xmin + i*desc.xstep,
			       desc.ymin + j*desc.ystep,
			       desc.zmin + k*desc.zstep );
	    }
	}
    }
    // Elements
    HexElem elem;
    for(int k=0;k<nelemz;++k) {
	for(int j=0;j<nelemz;++j) {
	    for(int i=0;i<nelemz;++i) {
		elem.v[0] = pointidx(i  ,j  ,k  ,nelemx,nelemy,nelemz);
		elem.v[1] = pointidx(i+1,j  ,k  ,nelemx,nelemy,nelemz);
		elem.v[2] = pointidx(i+1,j+1,k  ,nelemx,nelemy,nelemz);
		elem.v[3] = pointidx(i  ,j+1,k  ,nelemx,nelemy,nelemz);
		elem.v[4] = pointidx(i  ,j  ,k+1,nelemx,nelemy,nelemz);
		elem.v[5] = pointidx(i+1,j  ,k+1,nelemx,nelemy,nelemz);
		elem.v[6] = pointidx(i+1,j+1,k+1,nelemx,nelemy,nelemz);
		elem.v[7] = pointidx(i  ,j+1,k+1,nelemx,nelemy,nelemz);
		mesh.add_elem(0, elem);
	    }
	}
    }
}


void output_mesh_part_txt(MeshPart& loc, Mesh3D& mesh)
{
    int i;
    int part = loc.part;

    printf("%04d : number of elements = %d\n", part, loc.n_elems());
    char fname[2048];
    snprintf(fname, sizeof(fname), "mesh4spec.%04d", part);
    FILE* of = fopen(fname, "w");

    fprintf(of,"     3\n"); // Dims
    fprintf(of," F\n");     // Solid fluid
    fprintf(of," F\n");     // All fluid
    fprintf(of," F\n");     // Neumann
    fprintf(of," Local nodes\n");
    fprintf(of,"%6d\n", loc.n_nodes());
    fprintf(of," F\n");     // Curve
    for(int k=0; k<loc.n_nodes();++k) {
	int g = loc.n_l2g[k];
	fprintf(of, " %20f      %20f      %20f\n", mesh.m_xco[g], mesh.m_yco[g], mesh.m_zco[g] );
    }
    fprintf(of," Nb elements\n");
    fprintf(of,"%6d\n", loc.n_elems());
    fprintf(of," Materials\n");
    fprintf(of,"%6d\n", mesh.n_materials());
    for(int k=0; k<loc.n_elems(); ++k) {
	int g = loc.e_l2g[k];
	fprintf(of, "%6d       F\n", mesh.m_mat[g]);
    }
    fprintf(of," Global nodes for elements\n");    fprintf(of,"%6d\n", mesh.nodes_per_elem() );
    for(int k=0;k<loc.n_elems();++k) {
	int g = loc.e_l2g[k];
	for(i=mesh.m_elems_offs[g];i<mesh.m_elems_offs[g+1];++i) {
	    int ln = loc.n_g2l[mesh.m_elems[i]];
	    fprintf(of, "%6d", ln);
	}
	fprintf(of, "\n");
    }
    fprintf(of, "Faces\n");
    fprintf(of, "%6d\n", loc.n_faces());
    for(int k=0;k<loc.n_elems();++k) {
	fprintf(of, "%6d%6d%6d%6d%6d%6d\n",
		loc.faces[6*k+0],
		loc.faces[6*k+1],
		loc.faces[6*k+2],
		loc.faces[6*k+3],
		loc.faces[6*k+4],
		loc.faces[6*k+5]);
	fprintf(of, "%6d%6d%6d%6d%6d%6d\n",
		loc.faces_orient[6*k+0],
		loc.faces_orient[6*k+1],
		loc.faces_orient[6*k+2],
		loc.faces_orient[6*k+3],
		loc.faces_orient[6*k+4],
		loc.faces_orient[6*k+5]);
    }
}

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


int main(int argc, char**argv)
{
    Mesh3D mesh;
    RectMesh desc;

    desc.xmin = -100;
    desc.xmax =  100;
    desc.ymin = -100;
    desc.ymax =  100;
    desc.zmin = -100;
    desc.zmax =  100;
    desc.xstep = desc.ystep = desc.zstep = 10;
    desc.pmls.N = false;
    desc.pmls.S = false;
    desc.pmls.E = false;
    desc.pmls.W = false;
    desc.pmls.U = false;
    desc.pmls.D = true;

    mesh.add_material();

    init_rectangular_mesh(desc, mesh);
    int NPROCS=2;
    mesh.partition_mesh(NPROCS);

    for(int part=0;part<NPROCS;++part) {
	MeshPart loc;
	mesh.compute_local_part(part, loc);
	output_mesh_part_h5(loc, mesh);
	output_mesh_part_xmf(loc, mesh);
    }
    output_all_meshes_xmf(NPROCS);
    for(int part=0;part<NPROCS;++part) {
	output_mesh_part_h5_comm(part, mesh);
    }
}
// Local Variables:
// mode: c++
// c-file-style: "stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=4 ts=8 tw=80 smartindent */
