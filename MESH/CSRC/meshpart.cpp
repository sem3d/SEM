
#include <cstdio>
#include <cassert>
#include <vector>
#include "mesh.h"
#include "meshpart.h"

using std::vector;


/// Reference numbering for faces
static FaceDesc RefFace[6] = {
    { {0,1,2,3}, {0, 1, 2, 3} },
    { {0,1,5,4}, {0, 4, 5, 6} },
    { {1,2,6,5}, {1, 7, 8, 4} },
    { {3,2,6,7}, {2, 7, 9,10} },
    { {0,3,7,4}, {3,10,11, 6} },
    { {4,5,6,7}, {5, 8, 9,11} },
};

/// Reference numbering of edges
static int RefEdge[12][2] = {
    { 0, 1},
    { 1, 2},
    { 3, 2},
    { 0, 3},
    { 1, 5},
    { 4, 5},
    { 0, 4},
    { 2, 6},
    { 5, 6},
    { 7, 6},
    { 3, 7},
    { 4, 7}
};

void Mesh3DPart::compute_part()
{
    /// Handle all elements on this node and those that touches it
    for(int k=0;k<m_mesh.n_elems();++k) {
        if (m_mesh.elem_part(k)==m_proc) {
            handle_local_element(k);
        } else {
            handle_neighbour_element(k);
        }
    }
    printf("Created %ld facets\n", m_face_to_id.size());
}

int Mesh3DPart::add_facet(int n[4], int dom)
{
    int nf;
    face_map_t::iterator it;
    PFace facet(n, dom);
    it = m_face_to_id.find(facet);
    if (it==m_face_to_id.end()) {
        // New face
        nf = m_face_to_id.size();
        m_face_to_id[facet] = nf;
    } else {
        nf = it->second;
    }
    return nf;
}

int Mesh3DPart::add_edge(int v0, int v1, int dom)
{
    int ne;
    edge_map_t::iterator it;
    PEdge edge(v0, v1, dom);
    it = m_edge_to_id.find(edge);
    if (it==m_edge_to_id.end()) {
        // New edge
        ne = m_edge_to_id.size();
        m_edge_to_id[edge] = ne;
    } else {
        ne = it->second;
    }
    return ne;
}

int Mesh3DPart::add_vertex(int v0, int dom)
{
    int nv;
    vertex_map_t::iterator it;
    PVertex vertex(v0, dom);
    it = m_vertex_to_id.find(vertex);
    if (it==m_vertex_to_id.end()) {
        // New vertex
        nv = m_vertex_to_id.size();
        m_vertex_to_id[vertex] = nv;
    } else {
        nv = it->second;
    }
    return nv;
}

int Mesh3DPart::add_node(int v0)
{
    int nv;
    std::map<int,int>::iterator it;
    it = m_nodes_to_id.find(v0);
    if (it==m_nodes_to_id.end()) {
        // New vertex
        nv = m_nodes_to_id.size();
        m_nodes_to_id[v0] = nv;
    } else {
        nv = it->second;
    }
    return nv;
}

void Mesh3DPart::handle_local_element(int el)
{
    int e0 = m_mesh.m_elems_offs[el];
    int dom = m_mesh.get_elem_domain(el);

    m_elems.push_back(el);
    // Assign all 6 faces
    for(int fc=0;fc<6;++fc) {
        int n[4];
        for(int p=0;p<4;++p) {
            n[p] = m_mesh.m_elems[e0 + RefFace[fc].v[p]];
        }
        int nf = add_facet(n, dom);
        m_elems_faces.push_back(nf);
    }
    for(int ed=0;ed<12;++ed) {
        int v0 = RefEdge[ed][0];
        int v1 = RefEdge[ed][1];
        int ne = add_edge(m_mesh.m_elems[e0 + v0], m_mesh.m_elems[e0 + v1], dom);
        m_elems_edges.push_back(ne);
    }
    for(int vx=0;vx<8;++vx) {
        int gid = m_mesh.m_elems[e0 + vx];
        int vid = add_vertex(gid, dom);
        add_node(gid);
        m_elems_vertices.push_back(vid);
    }
    for(int vx=8;vx<m_mesh.nodes_per_elem();++vx) {
        add_node(m_mesh.m_elems[e0 + vx]);
    }
}

void Mesh3DPart::handle_neighbour_element(int el)
{
    // use adjacency map to tell if a neigbouring element belong to this proc
    // this works only if we used ncommon=1 in MeshToDual
    bool found = false;
    int contact;
    for(int k=m_mesh.m_xadj[el];k<m_mesh.m_xadj[el+1];++k) {
        int neighbour = m_mesh.m_adjncy[k];
        if (m_mesh.elem_part(neighbour)==m_proc) {
            found = true;
            break;
        }
    }
    if (!found) return;

    int source_proc = m_mesh.elem_part(el);
    // We have a neighbouring element
    // add only those facets, edge, vertice that are in common with our proc
    int e0 = m_mesh.m_elems_offs[el];
    int dom = m_mesh.get_elem_domain(el);
    int share_pt[8] = {0, 0, 0, 0, 0, 0, 0, 0 };
    vector<int> elems;;
    for(int k=0;k<8;++k) {
        int vertex_id = m_mesh.m_elems[e0+k];
        elems.clear();
        m_mesh.m_vertex_to_elem.vertex_to_elements(vertex_id, elems);
        for(size_t n=0;n<elems.size();++n) {
            if (m_mesh.elem_part(elems[n])==m_proc) {
                share_pt[k] = 1;
                break;
            }
        }
    }
    // Assign all 6 faces
    for(int fc=0;fc<6;++fc) {
        int n[4];
        contact = 0;
        for(int p=0;p<4;++p) {
            int vx = RefFace[fc].v[p];
            n[p] = m_mesh.m_elems[e0 + vx];
            contact += share_pt[vx];
        }
        if (contact==4) {
            int nf = add_facet(n, dom);
            m_comm[source_proc].m_faces.push_back(nf);
        } else if (contact!=0) {
            //printf("!!Face from part %d with %d nodes in contact with part %d\n", source_proc, contact, m_proc);
        }
    }
    for(int ed=0;ed<12;++ed) {
        int v0 = RefEdge[ed][0];
        int v1 = RefEdge[ed][1];
        contact = share_pt[v0] + share_pt[v1];
        if (contact==2) {
            int ne = add_edge(m_mesh.m_elems[e0 + v0], m_mesh.m_elems[e0 + v1], dom);
            m_comm[source_proc].m_edges.push_back(ne);
        } else if (contact!=0) {
            //printf("!!Edge from part %d with %d nodes in contact with part %d\n", source_proc, contact, m_proc);
        }
    }
    for(int vx=0;vx<8;++vx) {
        if (share_pt[vx]==1) {
            int nv = add_vertex(m_mesh.m_elems[e0 + vx], dom);
            m_comm[source_proc].m_vertices.push_back(nv);
        }
    }
}

void Mesh3DPart::get_local_nodes(std::vector<double>& nodes) const
{
    std::map<int,int>::const_iterator it;
    nodes.resize(3*m_nodes_to_id.size());
    for(it=m_nodes_to_id.begin();it!=m_nodes_to_id.end();++it) {
        int gid = it->first;
        int lid = it->second;
        nodes[3*lid+0] = m_mesh.m_xco[gid];
        nodes[3*lid+1] = m_mesh.m_yco[gid];
        nodes[3*lid+2] = m_mesh.m_zco[gid];
    }
}

void Mesh3DPart::get_local_elements(std::vector<int>& elems) const
{
    std::map<int,int>::const_iterator it;
    elems.resize(m_elems.size()*8);

    for(size_t k=0;k<m_elems.size();++k) {
        int el = m_elems[k];
        int e0 = m_mesh.m_elems_offs[el];
        for(int n=0;n<m_mesh.nodes_per_elem();++n) {
            int gid = m_mesh.m_elems[e0+n];
            int lid = get(m_nodes_to_id, gid, -1);
            elems[8*k + n] = lid;
        }
    }
}

void Mesh3DPart::get_local_materials(std::vector<int>& mats) const
{
    mats.resize(m_elems.size()*2);

    for(size_t k=0;k<m_elems.size();++k) {
        int el = m_elems[k];
        mats[2*k+0] = m_mesh.m_mat[el];
        mats[2*k+1] = 0;
    }

}

void Mesh3DPart::get_local_faces(std::vector<int>& faces) const
{
    faces.resize(m_face_to_id.size()*4);
    face_map_t::const_iterator it;
    int fc=0;
    for(it=m_face_to_id.begin();it!=m_face_to_id.end();++it) {
        fc = it->second;
        for(int p=0;p<4;++p) {
            int gid = it->first.n[p];
            int lid = get(m_nodes_to_id, gid, -1);
            faces[4*fc + p] = lid;
        }
    }
}

void Mesh3DPart::get_local_edges(std::vector<int>& edges) const
{
    edges.resize(m_edge_to_id.size()*2);
    edge_map_t::const_iterator it;
    int ed=0;
    for(it=m_edge_to_id.begin();it!=m_edge_to_id.end();++it) {
        ed = it->second;
        for(int p=0;p<2;++p) {
            int gid = it->first.n[p];
            int lid = get(m_nodes_to_id, gid, -1);
            edges[2*ed + p] = lid;
        }
    }
}

void Mesh3DPart::output_mesh_part()
{
    char fname[2048];
    vector<double> tmpd;
    vector<int> tmpi, tmpi1;

    printf("%04d : number of elements = %d\n", m_proc, n_elems());
    snprintf(fname, sizeof(fname), "mesh4spec.%04d.h5", m_proc);
    hid_t fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    h5h_create_attr(fid, "ndims", 3);
    h5h_create_attr(fid, "n_processors", m_mesh.n_parts());
    h5h_create_attr(fid, "n_materials", m_mesh.n_materials());
    h5h_create_attr(fid, "n_elements", n_elems());
    h5h_create_attr(fid, "n_faces", n_faces());
    h5h_create_attr(fid, "n_edges", n_edges());
    h5h_create_attr(fid, "n_vertices", n_vertices());
    h5h_create_attr(fid, "solid_fluid", false);
    h5h_create_attr(fid, "solid_fluid_loc", false);
    h5h_create_attr(fid, "all_fluid", false);
    h5h_create_attr(fid, "neumann_present", false);
    h5h_create_attr(fid, "neumann_present_loc", false);
    h5h_create_attr(fid, "curve", false);

    get_local_nodes(tmpd);
    h5h_write_dset_2d(fid, "local_nodes", n_nodes(), 3, &tmpd[0]);

    get_local_elements(tmpi);
    h5h_write_dset_2d(fid, "elements", n_elems(), 8, &tmpi[0]);

    get_local_materials(tmpi);
    h5h_write_dset_2d(fid, "material", n_elems(), 2, &tmpi[0]);
    //
    h5h_write_dset_2d(fid, "faces", n_elems(), 6, &m_elems_faces[0]);
    get_local_faces(tmpi);
    h5h_write_dset_2d(fid, "faces_def", n_elems(), 4, &tmpi[0]);
    //
    h5h_write_dset_2d(fid, "edges", n_elems(), 12, &m_elems_edges[0]);
    get_local_edges(tmpi);
    h5h_write_dset_2d(fid, "edges_def", n_elems(), 2, &tmpi[0]);
    //
    h5h_write_dset_2d(fid, "vertices", n_elems(), 8, &m_elems_vertices[0]);
    //h5h_write_dset(fid, "vertices_to_global", n_l2g.size(), &n_l2g[0]);

    //
    H5Fclose(fid);
}

void Mesh3DPart::output_mesh_part_xmf()
{
    char fname[2048];
    FILE* f;

    snprintf(fname, sizeof(fname), "mesh4spec.%04d.xmf", m_proc);
    f = fopen(fname,"w");
    fprintf(f, "<?xml version=\"1.0\" ?>\n");
    fprintf(f, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n");
    fprintf(f, "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n");
    fprintf(f, "  <Domain>\n");
    fprintf(f, "    <Grid name=\"mesh.%04d\">\n", m_proc);
    fprintf(f, "      <Topology Type=\"Hexahedron\" NumberOfElements=\"%d\">\n", n_elems());
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%d 8\">\n", n_elems());
    fprintf(f, "mesh4spec.%04d.h5:/elements\n",m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Topology>\n");
    fprintf(f, "      <Geometry Type=\"XYZ\">\n");
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%d 3\">\n", n_nodes());
    fprintf(f, "mesh4spec.%04d.h5:/local_nodes\n", m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Geometry>\n");
    fprintf(f, "    </Grid>\n");
/*
    fprintf(f, "<Attribute Name=\"Mat\" Center=\"Cell\" AttributeType=\"Vector\" Dimensions=\"    44585 3\">");
    fprintf(f, "<DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"    44585 3\">");
    fprintf(f, "Rsem0001/sem_field.0000.h5:/displ");
    fprintf(f, "</DataItem>");
    fprintf(f, "</Attribute>");

*/
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
