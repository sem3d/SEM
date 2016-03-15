
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include "mesh.h"
#include "meshpart.h"

using std::vector;



void Mesh3DPart::compute_part()
{
    /// Handle all elements on this node and those that touches it
    for(int k=0;k<m_mesh.n_elems();++k) {
        if (m_mesh.elem_part(k)==m_proc) {
            bool border = is_border_element(k);
            handle_local_element(k, border);
        }
    }
    compute_face_communications();
    compute_edge_communications();
    compute_vertex_communications();
    // Walk the surfaces and keep what belongs to us...
    std::map<std::string, Surface*>::const_iterator it;
    for(it=m_mesh.m_surfaces.begin();it!=m_mesh.m_surfaces.end();++it) {
        handle_surface(it->second);
    }
    printf("Created %ld facets\n", m_face_to_id.size());
}


bool Mesh3DPart::is_border_element(int el)
{
    for(int k=m_mesh.m_xadj[el];k<m_mesh.m_xadj[el+1];++k) {
        int neighbour = m_mesh.m_adjncy[k];
        if (m_mesh.elem_part(neighbour)!=m_proc) {
            return true;
        }
    }
    return false;
}

void Mesh3DPart::get_neighbour_elements(int nn, const int* n, std::set<int>& elemset)
{
    std::set<int> elems, temp;
    m_mesh.m_vertex_to_elem.vertex_to_elements(n[0], elemset);
    for(int k=1;k<nn;++k) {
        int vertex_id = n[k];
        elems.clear();
        temp.clear();
        m_mesh.m_vertex_to_elem.vertex_to_elements(n[k], elems);
        std::set_intersection(elemset.begin(), elemset.end(),
                              elems.begin(), elems.end(),
                              std::inserter(temp, temp.begin()));
        elemset.swap(temp);
    }
}

void Mesh3DPart::compute_face_communications()
{
    // Move through all border vertex,edges, faces
    // and add them as communications if they match any
    // neighbour's from other processors
    // if two edge (faces, vertex) from different domains
    // lie on different processor, then each proc must have
    // the two domains.
    face_map_t::const_iterator it;
    // gross, but makes sure that we iterate over a stable container
    face_map_t face_to_id(m_face_to_id);
    for(it=face_to_id.begin();it!=face_to_id.end();++it) {
        if (!m_face_border[it->second]) continue;
        const PFace& fc = it->first;
        int nfc = it->second;
        std::set<int> elemset;
        get_neighbour_elements(4, fc.n, elemset);
        std::set<int>::const_iterator it;
        for(it=elemset.begin();it!=elemset.end();++it) {
            if (m_mesh.elem_part(*it)==m_proc) continue;
            handle_neighbouring_face(nfc, fc, *it);
        }
    }
}

void Mesh3DPart::compute_edge_communications()
{
    edge_map_t::const_iterator it;
    edge_map_t edge_to_id(m_edge_to_id);
    for(it=edge_to_id.begin();it!=edge_to_id.end();++it) {
        if (!m_edge_border[it->second]) continue;
        const PEdge& ed = it->first;
        int ned = it->second;
        std::set<int> elemset;
        get_neighbour_elements(2, ed.n, elemset);
        std::set<int>::const_iterator it;
        for(it=elemset.begin();it!=elemset.end();++it) {
            if (m_mesh.elem_part(*it)==m_proc) continue;
            handle_neighbouring_edge(ned, ed, *it);
        }
    }
}

void Mesh3DPart::compute_vertex_communications()
{
    vertex_map_t::const_iterator it;
    vertex_map_t vertex_to_id(m_vertex_to_id);
    for(it=vertex_to_id.begin();it!=vertex_to_id.end();++it) {
        if (!m_vertex_border[it->second]) continue;
        const PVertex& vx = it->first;
        int nvx = it->second;
        std::set<int> elemset;
        get_neighbour_elements(1, vx.n, elemset);
        std::set<int>::const_iterator it;
        for(it=elemset.begin();it!=elemset.end();++it) {
            if (m_mesh.elem_part(*it)==m_proc) continue;
            handle_neighbouring_vertex(nvx, vx, *it);
        }
    }
}
void Mesh3DPart::handle_neighbouring_face(int lnf, const PFace& fc, int el)
{
    int dom = m_mesh.get_elem_domain(el);
    int e0 = m_mesh.m_elems_offs[el];
    const int *eldef = &m_mesh.m_elems[e0];
    int source_proc = m_mesh.elem_part(el);
    for(int nf=0;nf<6;++nf) {
        int n[4];
        for(int i=0;i<4;++i) n[i] = eldef[RefFace[nf].v[i]];
        PFace nfc(n,dom);
        if (!nfc.eq_geom(fc)) continue;
        /// add local facet
        m_comm[source_proc].m_faces[fc]=lnf;
        if (dom!=fc.domain()) {
            // if neighbour's domain is different add a new facet
            int inf = add_facet(nfc, true);
            m_comm[source_proc].m_faces[nfc]=inf;
        }
    }
}

void Mesh3DPart::handle_neighbouring_edge(int lne, const PEdge& ed, int el)
{
    int e0 = m_mesh.m_elems_offs[el];
    const int *eldef = &m_mesh.m_elems[e0];
    int source_proc = m_mesh.elem_part(el);
    for(int ne=0;ne<12;++ne) {
        int n[2];
        for(int i=0;i<2;++i) n[i] = eldef[RefEdge[ne][i]];
        PEdge ned(n[0], n[1],0);
        if (!ned.eq_geom(ed)) continue;
        // We are guaranteed that the edge exists locally and on the remote
        // We add all domains common to the two vertices of the edge since
        // a domain not local to the processor and the source_proc may still
        // be created on source_proc by a third processor

        int dommask = m_mesh.m_vertex_domains[n[0]] & m_mesh.m_vertex_domains[n[1]];
        for(int dom=0;dom<=DM_MAX;++dom) {
            if ((dommask & (1<<dom))==0) continue;
            PEdge loced(n[0], n[1], dom);
            int ne = add_edge(loced, true);
            m_comm[source_proc].m_edges[loced] = ne;
        }
    }
}

void Mesh3DPart::handle_neighbouring_vertex(int lnv, const PVertex& vx, int el)
{
    int e0 = m_mesh.m_elems_offs[el];
    const int *eldef = &m_mesh.m_elems[e0];
    int source_proc = m_mesh.elem_part(el);
    int dommask = m_mesh.m_vertex_domains[vx.n[0]];
    for(int dom=0;dom<=DM_MAX;++dom) {
        if ((dommask & (1<<dom))==0) continue;
        PVertex locvx(vx.n[0], dom);
        int nv = add_vertex(locvx, true);
        m_comm[source_proc].m_vertices[locvx] = nv;
    }
}

void Mesh3DPart::handle_surface(const Surface* surf)
{
    // brutal, test each face if present...
    Surface* new_surf = new Surface(surf->name());
    face_map_t::const_iterator itfound, it;
    for(it=surf->m_faces.begin();it!=surf->m_faces.end();++it) {
        itfound = m_face_to_id.find(it->first);
        if (itfound!=m_face_to_id.end()) {
            const PFace& fc = it->first;
            new_surf->add_face(fc, itfound->second);
            int dom = fc.domain();
            for(int k=0;k<4;++k) {
                int va = fc.n[k];
                int vb = fc.n[(k+1)%4];
                PEdge ed(va, vb, dom);
                int eid = get(m_edge_to_id, ed, -1);
                assert(eid!=-1);
                new_surf->add_edge(ed, eid);
                PVertex vx(va, dom);
                int vid = get(m_vertex_to_id, vx, -1);
                assert(vid!=-1);
                new_surf->add_vertex(vx, vid);
            }
        }
    }
    // We keep all surfaces even if locally empty
    m_surfaces.push_back(new_surf);
}

int Mesh3DPart::add_facet(const PFace& facet, bool border)
{
    int nf;
    face_map_t::iterator it;
    it = m_face_to_id.find(facet);
    if (it==m_face_to_id.end()) {
        // New face
        nf = m_face_to_id.size();
        m_face_border.push_back(border);
        m_face_to_id[facet] = nf;
    } else {
        nf = it->second;
        m_face_border[nf] = border || m_face_border[nf];
    }
    return nf;
}

int Mesh3DPart::add_edge(const PEdge& edge, bool border)
{
    int ne;
    edge_map_t::iterator it;
    it = m_edge_to_id.find(edge);
    if (it==m_edge_to_id.end()) {
        // New edge
        ne = m_edge_to_id.size();
        m_edge_border.push_back(border);
        m_edge_to_id[edge] = ne;
    } else {
        ne = it->second;
        m_edge_border[ne] = border || m_edge_border[ne];
    }
    return ne;
}

int Mesh3DPart::add_vertex(const PVertex& vertex, bool border)
{
    int nv;
    vertex_map_t::iterator it;
    it = m_vertex_to_id.find(vertex);
    if (it==m_vertex_to_id.end()) {
        // New vertex
        nv = m_vertex_to_id.size();
        m_vertex_border.push_back(border);
        m_vertex_to_id[vertex] = nv;
    } else {
        nv = it->second;
        m_vertex_border[nv] = border || m_vertex_border[nv];
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

void Mesh3DPart::handle_local_element(int el, bool is_border)
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
        PFace facet(n, dom);
        int nf = add_facet(facet, is_border);
        m_elems_faces.push_back(nf);
    }
    for(int ed=0;ed<12;++ed) {
        int v0 = RefEdge[ed][0];
        int v1 = RefEdge[ed][1];
        PEdge edge(m_mesh.m_elems[e0 + v0], m_mesh.m_elems[e0 + v1], dom);
        int ne = add_edge(edge, is_border);
        m_elems_edges.push_back(ne);
    }
    for(int vx=0;vx<8;++vx) {
        int gid = m_mesh.m_elems[e0 + vx];
        PVertex vertex(gid, dom);
        int vid = add_vertex(vertex, is_border);
        add_node(gid);
        m_elems_vertices.push_back(vid);
    }
    for(int vx=8;vx<m_mesh.nodes_per_elem();++vx) {
        int gid = m_mesh.m_elems[e0 + vx];
        add_node(gid);
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
    int nctl_nodes = m_mesh.nodes_per_elem();
    elems.resize(m_elems.size()*nctl_nodes);

    for(size_t k=0;k<m_elems.size();++k) {
        int el = m_elems[k];
        int e0 = m_mesh.m_elems_offs[el];
        for(int n=0;n<nctl_nodes;++n) {
            int gid = m_mesh.m_elems[e0+n];
            int lid = get(m_nodes_to_id, gid, -1);
            elems[nctl_nodes*k + n] = lid;
        }
    }
}

void Mesh3DPart::get_local_materials(std::vector<int>& mats, std::vector<int>& doms) const
{
    mats.resize(m_elems.size());
    doms.resize(m_elems.size());

    for(size_t k=0;k<m_elems.size();++k) {
        int el = m_elems[k];
        mats[k] = m_mesh.m_mat[el];
        int dom = m_mesh.get_elem_domain(el);
        doms[k] = dom;
        assert(dom>0 && dom<=DM_MAX);
    }

}

void Mesh3DPart::get_local_faces(std::vector<int>& faces, std::vector<int>& doms) const
{
    faces.resize(m_face_to_id.size()*4);
    doms.resize(m_face_to_id.size());
    face_map_t::const_iterator it;
    int fc=0;
    for(it=m_face_to_id.begin();it!=m_face_to_id.end();++it) {
        fc = it->second;
        for(int p=0;p<4;++p) {
            int gid = it->first.n[p];
            int lid = get(m_nodes_to_id, gid, -1);
            faces[4*fc + p] = lid;
            doms[fc] = it->first.domain();
        }
    }
}

void Mesh3DPart::get_local_edges(std::vector<int>& edges, std::vector<int>& doms) const
{
    edges.resize(2*m_edge_to_id.size());
    doms.resize(m_edge_to_id.size());
    edge_map_t::const_iterator it;
    int ed=0;
    for(it=m_edge_to_id.begin();it!=m_edge_to_id.end();++it) {
        ed = it->second;
        for(int p=0;p<2;++p) {
            int gid = it->first.n[p];
            int lid = get(m_nodes_to_id, gid, -1);
            edges[2*ed + p] = lid;
            doms[ed] = it->first.domain();
        }
    }
}

void Mesh3DPart::get_local_vertices_dom(std::vector<int>& doms) const
{
    doms.resize(m_vertex_to_id.size());
    vertex_map_t::const_iterator it;
    int vx=0;
    for(it=m_vertex_to_id.begin();it!=m_vertex_to_id.end();++it) {
        vx = it->second;
        for(int p=0;p<2;++p) {
            doms[vx] = it->first.domain();
        }
    }
}

void Mesh3DPart::global_to_local_ids(std::vector<int>& ids) const
{
    std::map<int,int>::const_iterator it, itend;
    itend = m_nodes_to_id.end();
    for(int k=0;k<ids.size();++k) {
        it = m_nodes_to_id.find(ids[k]);
        if (it!=itend)
            ids[k] = it->second;
    }
}

void Mesh3DPart::get_face_coupling(int d0, int d1, std::vector<int>& cpl, std::vector<int>& orient) const
{
    bool swapped=false;
    cpl.clear();
    orient.clear();
    assert(d0!=d1);
    if (d0>d1) {
        swapped = true; // preserve the domains in the required order
        std::swap(d0,d1);
    }
    face_map_t::const_iterator it0, it;

    it=m_face_to_id.begin();
    it0=m_face_to_id.begin();
    if (it==m_face_to_id.end()) return;
    ++it;
    if (it==m_face_to_id.end()) return;
    // it0 and it are two consecutive faces
    for(;it!=m_face_to_id.end();++it,++it0) {
        if (!it->first.eq_geom(it0->first)) continue;
        // We have two faces equal with different domain, 
        if (it0->first.domain()!=d0) continue;
        if (it->first.domain()!=d1) continue;
        if (swapped) {
            cpl.push_back(it->second);
            cpl.push_back(it0->second);
            orient.push_back(it->first.orient);
            orient.push_back(it0->first.orient);
        } else {
            cpl.push_back(it0->second);
            cpl.push_back(it->second);
            orient.push_back(it0->first.orient);
            orient.push_back(it->first.orient);
        }
    }
}

void Mesh3DPart::get_edge_coupling(int d0, int d1, std::vector<int>& cpl) const
{
    bool swapped=false;
    cpl.clear();
    assert(d0!=d1);
    if (d0>d1) {
        swapped = true;
        std::swap(d0,d1);
    }
    edge_map_t::const_iterator it0, it;

    for(it0=m_edge_to_id.begin();it0!=m_edge_to_id.end();++it0) {
        it = it0;
        ++it;
        if (it0->first.domain()!=d0) continue;
        while(it!=m_edge_to_id.end()) {
            if (!it->first.eq_geom(it0->first)) break;
            // We have two edges equal with different domain,
            if (it->first.domain()!=d1) {
                ++it;
                continue;
            }
            if (swapped) {
                cpl.push_back(it->second);
                cpl.push_back(it0->second);
            } else {
                cpl.push_back(it0->second);
                cpl.push_back(it->second);
            }
            ++it;
        }
    }
}

void Mesh3DPart::get_vertex_coupling(int d0, int d1, std::vector<int>& cpl) const
{
    bool swapped=false;
    cpl.clear();
    assert(d0!=d1);
    if (d0>d1) {
        swapped = true;
        std::swap(d0,d1);
    }
    vertex_map_t::const_iterator it0, it;

    for(it0=m_vertex_to_id.begin();it0!=m_vertex_to_id.end();++it0) {
        it = it0;
        ++it;
        if (it0->first.domain()!=d0) continue;
        while(it!=m_vertex_to_id.end()) {
            if (!it->first.eq_geom(it0->first)) break;
            // We have two vertexs equal with different domain,
            if (it->first.domain()!=d1) {
                ++it;
                continue;
            }
            if (swapped) {
                cpl.push_back(it->second);
                cpl.push_back(it0->second);
            } else {
                cpl.push_back(it0->second);
                cpl.push_back(it->second);
            }
            ++it;
        }
    }
}

void Mesh3DPart::output_mesh_attributes(hid_t fid)
{
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
}

void Mesh3DPart::output_local_mesh(hid_t fid)
{
    int elem_doms[5] = {0,0,0,0,0};
    int face_doms[5] = {0,0,0,0,0};
    int edge_doms[5] = {0,0,0,0,0};
    int vert_doms[5] = {0,0,0,0,0};
    vector<double> tmpd;
    vector<int> tmpi, tmpi1;

    get_local_nodes(tmpd);
    h5h_write_dset_2d(fid, "local_nodes", n_nodes(), 3, &tmpd[0]);

    get_local_elements(tmpi);
    h5h_write_dset_2d(fid, "elements", n_elems(), m_mesh.nodes_per_elem(), &tmpi[0]);

    get_local_materials(tmpi,tmpi1);
    h5h_write_dset(fid, "material", n_elems(), &tmpi[0]);
    h5h_write_dset(fid, "domains", n_elems(), &tmpi1[0]);
    for(unsigned k=0;k<tmpi1.size();++k) {
        elem_doms[tmpi1[k]]++;
    }
    //
    h5h_write_dset_2d(fid, "faces", n_elems(), 6, &m_elems_faces[0]);
    get_local_faces(tmpi, tmpi1);
    h5h_write_dset_2d(fid, "faces_def", n_faces(), 4, &tmpi[0]);
    h5h_write_dset(fid, "faces_dom", n_faces(), &tmpi1[0]);
    for(unsigned k=0;k<tmpi1.size();++k) {
        face_doms[tmpi1[k]]++;
    }
    //
    h5h_write_dset_2d(fid, "edges", n_elems(), 12, &m_elems_edges[0]);
    get_local_edges(tmpi, tmpi1);
    h5h_write_dset_2d(fid, "edges_def", n_edges(), 2, &tmpi[0]);
    h5h_write_dset(fid, "edges_dom", n_edges(), &tmpi1[0]);
    for(unsigned k=0;k<tmpi1.size();++k) {
        edge_doms[tmpi1[k]]++;
    }
    h5h_write_dset_2d(fid, "vertices", n_elems(), 8, &m_elems_vertices[0]);
    get_local_vertices_dom(tmpi);
    h5h_write_dset(fid, "vertices_dom", n_vertices(), &tmpi[0]);
    //
    printf("%04d : number of elements = (tot=%d/fpml=%d/spml=%d/fl=%d/sol=%d)\n", m_proc, n_elems(),
           elem_doms[1],elem_doms[2],elem_doms[3],elem_doms[4]);
    printf("%04d : number of faces    = (tot=%d/fpml=%d/spml=%d/fl=%d/sol=%d)\n", m_proc, n_faces(),
           face_doms[1],face_doms[2],face_doms[3],face_doms[4]);
    printf("%04d : number of edges    = (tot=%d/fpml=%d/spml=%d/fl=%d/sol=%d)\n", m_proc, n_edges(),
           edge_doms[1],edge_doms[2],edge_doms[3],edge_doms[4]);
}

void Mesh3DPart::write_surface_dom(hid_t gid, const Surface* surf, const char* pfx, int dom)
{
    vector<int> data, orient;
    char sface_data[100];
    char sface_orient[100];
    char sface_num [100];
    char sedge_data[100];
    char sedge_num [100];
    char svert_data[100];
    char svert_num [100];
    snprintf(sface_data, 100, "%s_faces", pfx);
    snprintf(sface_orient, 100, "%s_orient", pfx);
    snprintf(sface_num , 100, "n_%s_faces", pfx);
    snprintf(sedge_data, 100, "%s_edges", pfx);
    snprintf(sedge_num , 100, "n_%s_edges", pfx);
    snprintf(svert_data, 100, "%s_vertices", pfx);
    snprintf(svert_num , 100, "n_%s_vertices", pfx);

    surf->get_faces_data(dom, data, orient);
    h5h_create_attr(gid, sface_num, (int)data.size());
    h5h_write_dset(gid, sface_data, data);
    h5h_write_dset(gid, sface_orient, orient);

    surf->get_edges_data(dom, data);
    h5h_create_attr(gid, sedge_num, (int)data.size());
    h5h_write_dset(gid, sedge_data, data);

    surf->get_vertices_data(dom, data);
    h5h_write_dset(gid, svert_data, data);
    h5h_create_attr(gid, svert_num, (int)data.size());
}
void Mesh3DPart::output_surface(hid_t fid, const Surface* surf)
{
    std::vector<int> faces, orient;
    hid_t gid = H5Gcreate(fid, surf->name().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_surface_dom(gid, surf, "sl", DM_SOLID);
    write_surface_dom(gid, surf, "fl", DM_FLUID);
    write_surface_dom(gid, surf, "spml", DM_SOLID_PML);
    write_surface_dom(gid, surf, "fpml", DM_FLUID_PML);
    H5Gclose(gid);
}

void Mesh3DPart::output_mesh_part()
{
    char fname[2048];

    h5h_set_errhandler();

    printf("%04d : number of elements = %d\n", m_proc, n_elems());
    snprintf(fname, sizeof(fname), "mesh4spec.%04d.h5", m_proc);
    hid_t fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    output_mesh_attributes(fid);
    output_local_mesh(fid);

    hid_t surf_grp = H5Gcreate(fid, "Surfaces", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5h_create_attr(surf_grp, "n_surfaces", (int)m_surfaces.size());
    for(unsigned k=0;k<m_surfaces.size();++k) {
        output_surface(surf_grp, m_surfaces[k]);
    }
    H5Gclose(surf_grp);

    // Now write out inter-domain coupling
    // Couplage Solide-fluide :
    write_coupling_interface(fid, "sf", DM_SOLID, DM_FLUID);
    write_coupling_interface(fid, "sfpml", DM_SOLID_PML, DM_FLUID_PML);
    write_coupling_interface(fid, "fpml", DM_FLUID, DM_FLUID_PML);
    write_coupling_interface(fid, "spml", DM_SOLID, DM_SOLID_PML);

    // Write processors communications
    h5h_create_attr(fid, "tot_comm_proc", (int)m_comm.size());
    std::map<int,MeshPartComm>::const_iterator it;
    int k=0;
    for(it=m_comm.begin();it!=m_comm.end();++it,++k) {
        char gproc[200];
        snprintf(gproc, 200, "Proc%04d", k);
        hid_t gid = H5Gcreate(fid, gproc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        output_comm(gid, it->second, it->first);
    }

    H5Fclose(fid);
}


void Mesh3DPart::output_comm(hid_t gid, const MeshPartComm& comm, int dest)
{
    std::vector<int> tdefs, tids, tdoms;
    h5h_create_attr(gid, "proc_dest", dest);
    h5h_create_attr(gid, "n_faces", comm.n_faces());
    h5h_create_attr(gid, "n_edges", comm.n_edges());
    h5h_create_attr(gid, "n_vertices", comm.n_vertices());
    comm.get_faces(tdefs, tids, tdoms);
    h5h_write_dset(gid, "faces", tids);
    if (m_mesh.debug) {
        global_to_local_ids(tdefs);
        h5h_write_dset_2d(gid, "faces_def", 4, tdefs); // DEBUG
        h5h_write_dset(gid, "faces_dom", tdoms);
    }
    comm.get_edges(tdefs, tids, tdoms);
    h5h_write_dset(gid, "edges", tids);
    if (m_mesh.debug) {
        global_to_local_ids(tdefs);
        h5h_write_dset_2d(gid, "edges_def", 2, tdefs); // DEBUG
        h5h_write_dset(gid, "edges_dom", tdoms);
    }
    comm.get_vertices(tdefs, tids, tdoms);
    h5h_write_dset(gid, "vertices", tids);
    if (m_mesh.debug) {
        global_to_local_ids(tdefs);
        h5h_write_dset(gid, "vertices_defs", tdefs);
    }
    H5Gclose(gid);
    if (m_mesh.debug) {
        printf("Comm:%d->%d : F/E/V : (%d,%d,%d)\n", m_proc, dest,
               (int)comm.m_faces.size(), (int)comm.m_edges.size(),
               (int)comm.m_vertices.size());
    }

}
void Mesh3DPart::reorder_comm(MeshPartComm& comm)
{
}

void Mesh3DPart::write_coupling_interface(hid_t fid, const char* pfx, int d0, int d1)
{
    vector<int> tmpi, tmpo;
    char sface_data[100];
    char sface_orient[100];
    char sface_num [100];
    char sedge_data[100];
    char sedge_num [100];
    char svert_data[100];
    char svert_num [100];
    snprintf(sface_data, 100, "%s_faces", pfx);
    snprintf(sface_orient, 100, "%s_orient", pfx);
    snprintf(sface_num , 100, "n_%s_faces", pfx);
    snprintf(sedge_data, 100, "%s_edges", pfx);
    snprintf(sedge_num , 100, "n_%s_edges", pfx);
    snprintf(svert_data, 100, "%s_vertices", pfx);
    snprintf(svert_num , 100, "n_%s_vertices", pfx);

    get_face_coupling(d0, d1, tmpi, tmpo);
    h5h_create_attr(fid, sface_num, int(tmpi.size()/2) );
    h5h_write_dset_2d(fid, sface_data, tmpi.size()/2, 2, &tmpi[0]);
    h5h_write_dset_2d(fid, sface_orient, tmpo.size()/2, 2, &tmpo[0]);
    //
    get_edge_coupling(d0, d1, tmpi);
    h5h_create_attr(fid, sedge_num, int(tmpi.size()/2) );
    h5h_write_dset_2d(fid, sedge_data, tmpi.size()/2, 2, &tmpi[0]);
    //
    get_vertex_coupling(d0, d1, tmpi);
    h5h_create_attr(fid, svert_num, int(tmpi.size()/2) );
    h5h_write_dset_2d(fid, svert_data, tmpi.size()/2, 2, &tmpi[0]);
}

void Mesh3DPart::output_int_scalar(FILE* f, int indent, const char* aname, const char* atype, int n0, const char* field)
{
    char sind[50];
    for(int k=0;k<indent;++k) sind[k] = ' ';
    sind[indent] = 0;
    fprintf(f, "%s<Attribute Name=\"%s\" Center=\"%s\" AttributeType=\"Scalar\" Dimensions=\"%d\">\n", sind, aname, atype, n0);
    fprintf(f, "%s  <DataItem Format=\"HDF\" Datatype=\"Integer\" Dimensions=\"%d\">\n", sind, n0);
    fprintf(f, "mesh4spec.%04d.h5:%s\n", m_proc, field);
    fprintf(f, "%s  </DataItem>\n", sind);
    fprintf(f, "%s</Attribute>\n", sind);
}


void Mesh3DPart::output_int_constant(FILE* f, int indent, const char* aname, const char* atype, int val)
{
    char sind[50];
    for(int k=0;k<indent;++k) sind[k] = ' ';
    sind[indent] = 0;
    fprintf(f, "%s<Attribute Name=\"%s\" Center=\"%s\" AttributeType=\"Constant\" Dimensions=\"1\">\n", sind, aname, atype);
    fprintf(f, "%s  <DataItem Format=\"XML\" Datatype=\"Integer\" Dimensions=\"1\">\n", sind);
    fprintf(f, "%d\n", val);
    fprintf(f, "%s  </DataItem>\n", sind);
    fprintf(f, "%s</Attribute>\n", sind);
}

void Mesh3DPart::output_mesh_part_xmf()
{
    output_xmf_elements();
    output_xmf_faces();
    output_xmf_edges();
    output_xmf_vertices();
    output_xmf_comms();
}

void Mesh3DPart::output_xmf_header(FILE* f)
{
    fprintf(f, "<?xml version=\"1.0\" ?>\n");
    fprintf(f, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n");
    fprintf(f, "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n");
    fprintf(f, "  <Domain>\n");
}

void Mesh3DPart::output_xmf_footer(FILE* f)
{
    fprintf(f, "  </Domain>\n");
    fprintf(f, "</Xdmf>\n");
}

void Mesh3DPart::output_xmf_elements()
{
    char fname[2048];
    FILE* f;

    snprintf(fname, sizeof(fname), "mesh4spec.%04d.elems.xmf", m_proc);
    f = fopen(fname,"w");
    output_xmf_header(f);
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

    output_int_scalar(f, 6, "Mat", "Cell", n_elems(), "/material");
    output_int_scalar(f, 6, "Dom", "Cell", n_elems(), "/domains");
    fprintf(f, "    </Grid>\n");
    output_xmf_footer(f);
    fclose(f);

}


void Mesh3DPart::output_xmf_faces()
{
    char fname[2048];
    FILE* f;

    snprintf(fname, sizeof(fname), "mesh4spec.%04d.faces.xmf", m_proc);
    f = fopen(fname,"w");
    output_xmf_header(f);
    fprintf(f, "    <Grid name=\"faces.%04d\">\n", m_proc);
    fprintf(f, "      <Topology Type=\"Quadrilateral\" NumberOfElements=\"%d\">\n", n_faces());
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%d 4\">\n", n_faces());
    fprintf(f, "mesh4spec.%04d.h5:/faces_def\n",m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Topology>\n");
    fprintf(f, "      <Geometry Type=\"XYZ\">\n");
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%d 3\">\n", n_nodes());
    fprintf(f, "mesh4spec.%04d.h5:/local_nodes\n", m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Geometry>\n");
    output_int_scalar(f, 6, "FDom", "Cell", n_faces(), "/faces_dom");
    fprintf(f, "    </Grid>\n");
    output_xmf_footer(f);
    fclose(f);
}

void Mesh3DPart::output_xmf_edges()
{
    char fname[2048];
    FILE* f;

    snprintf(fname, sizeof(fname), "mesh4spec.%04d.edges.xmf", m_proc);
    f = fopen(fname,"w");
    output_xmf_header(f);
    fprintf(f, "    <Grid name=\"edges.%04d\">\n", m_proc);
    fprintf(f, "      <Topology Type=\"Polyline\" NumberOfElements=\"%d\">\n", n_edges());
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%d 2\">\n", n_edges());
    fprintf(f, "mesh4spec.%04d.h5:/edges_def\n",m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Topology>\n");
    fprintf(f, "      <Geometry Type=\"XYZ\">\n");
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%d 3\">\n", n_nodes());
    fprintf(f, "mesh4spec.%04d.h5:/local_nodes\n", m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Geometry>\n");
    output_int_scalar(f, 6, "FDom", "Cell", n_edges(), "/edges_dom");
    fprintf(f, "    </Grid>\n");
    output_xmf_footer(f);
    fclose(f);
}

void Mesh3DPart::output_xmf_vertices()
{
    char fname[2048];
    FILE* f;

    snprintf(fname, sizeof(fname), "mesh4spec.%04d.vertices.xmf", m_proc);
    f = fopen(fname,"w");
    output_xmf_header(f);

    output_xmf_footer(f);
    fclose(f);
}

void Mesh3DPart::output_xmf_comms()
{
    output_xmf_comms_faces();
    output_xmf_comms_edges();
}

void Mesh3DPart::output_xmf_comms_edges()
{
    char namebuf[2048];
    FILE* f;

    snprintf(namebuf, sizeof(namebuf), "mesh4spec.%04d.comms.edges.xmf", m_proc);
    f = fopen(namebuf,"w");
    output_xmf_header(f);

    fprintf(f, "  <Grid name=\"comms.%04d\" GridType=\"Collection\" CollectionType=\"Spatial\">\n", m_proc);
    int k=0;
    std::map<int,MeshPartComm>::const_iterator it;
    for(it=m_comm.begin();it!=m_comm.end();++it,++k) {
        int dest = it->first;

        // Edges
        fprintf(f, "    <Grid name=\"comms.%04d.%04d.e\">\n", m_proc, dest);
        int n_edges = it->second.n_edges();
        fprintf(f, "      <Topology Type=\"Polyline\" NumberOfElements=\"%d\">\n", n_edges);
        fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%d 2\">\n", n_edges);
        fprintf(f, "mesh4spec.%04d.h5:/Proc%04d/edges_def\n",m_proc, k);
        fprintf(f, "        </DataItem>\n");
        fprintf(f, "      </Topology>\n");
        fprintf(f, "      <Geometry Type=\"XYZ\">\n");
        fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%d 3\">\n", n_nodes());
        fprintf(f, "mesh4spec.%04d.h5:/local_nodes\n", m_proc);
        fprintf(f, "        </DataItem>\n");
        fprintf(f, "      </Geometry>\n");
        snprintf(namebuf, sizeof(namebuf), "/Proc%04d/edges_dom", k);
        output_int_scalar(f, 8, "Dom", "Cell", n_edges, namebuf);

        fprintf(f, "    </Grid>\n");
    }
    fprintf(f, "  </Grid>\n");
    output_xmf_footer(f);
    fclose(f);
}

void Mesh3DPart::output_xmf_comms_faces()
{
    char namebuf[2048];
    FILE* f;

    snprintf(namebuf, sizeof(namebuf), "mesh4spec.%04d.comms.faces.xmf", m_proc);
    f = fopen(namebuf,"w");
    output_xmf_header(f);

    fprintf(f, "  <Grid name=\"comms.%04d\" GridType=\"Collection\" CollectionType=\"Spatial\">\n", m_proc);
    int k=0;
    std::map<int,MeshPartComm>::const_iterator it;
    for(it=m_comm.begin();it!=m_comm.end();++it,++k) {
        int dest = it->first;

        // Faces
        fprintf(f, "    <Grid name=\"comms.%04d.%04d.e\">\n", m_proc, dest);
        int n_faces = it->second.n_faces();
        fprintf(f, "      <Topology Type=\"Quadrilateral\" NumberOfElements=\"%d\">\n", n_faces);
        fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%d 4\">\n", n_faces);
        fprintf(f, "mesh4spec.%04d.h5:/Proc%04d/faces_def\n",m_proc, k);
        fprintf(f, "        </DataItem>\n");
        fprintf(f, "      </Topology>\n");
        fprintf(f, "      <Geometry Type=\"XYZ\">\n");
        fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%d 3\">\n", n_nodes());
        fprintf(f, "mesh4spec.%04d.h5:/local_nodes\n", m_proc);
        fprintf(f, "        </DataItem>\n");
        fprintf(f, "      </Geometry>\n");
        snprintf(namebuf, sizeof(namebuf), "/Proc%04d/faces_dom", k);
        output_int_scalar(f, 8, "Dom", "Cell", n_faces, namebuf);
        fprintf(f, "    </Grid>\n");
    }
    fprintf(f, "  </Grid>\n");
    output_xmf_footer(f);
    fclose(f);
}





/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
