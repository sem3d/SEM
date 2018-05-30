
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include "mesh.h"
#include "meshpart.h"
#include "sem_gll.h"

using std::vector;

void Mesh3DPart::compute_part()
{
    /// Handle all elements on this node and those that touches it
    compute_gll();
    for(size_t k=0;k<m_mesh.n_elems();++k) {
        if (m_mesh.elem_part(k)==m_proc) {
            bool border = is_border_element(k);
            handle_local_element(k, border);
            //define_Eventual_Boundary_Surface( k );
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

void Mesh3DPart::define_Eventual_Boundary_Surface(int el)
{
#if 0
    int e0 = m_mesh.m_elems_offs[el];
    int dom = m_mesh.get_elem_domain(el);
    int nfc=-1;
    surf_info_map_t::const_iterator itt;
    itt=m_mesh.surfelem.find(el);
    if (itt!=m_mesh.surfelem.end()){
        int itmatname = itt->second.mat;
        int matdom    = itt->second.dom;
        Mesh3D & mesh=const_cast<Mesh3D&>(m_mesh);
        Surface* getsurface = mesh.get_surface(mesh.m_surf_matname[itmatname]);
        for (int fc=0; fc <6; ++fc) {
            bool add=true;
            int n[4];
            for(int p=0; p<4; ++p) { n[p] = m_mesh.m_elems[e0 + RefFace[fc].v[p]]; }
            for(int p=0; p<4; ++p) {
                if (std::find(itt->second.xxnodes.begin(),itt->second.xxnodes.end(), n[p]) ==
                    itt->second.xxnodes.end()) {
                    add=false;
                }
            }
            if (add) {
                nfc = fc; PFace face (n, dom);
                getsurface->add_face(face,m_mesh.m_materials[matdom].domain());
                break;
            }
            if ((!add)&&(fc>=5)&&(nfc==-1)) {
                printf("error: el %d, elf %d , e0 %d \n",el, itt->first, e0);
                exit(1);
            }
        }
        for (int ed=0; ed < 12; ++ed){
            int v0 = RefEdge[ed][0];
            int v1 = RefEdge[ed][1];
            int n[]={m_mesh.m_elems[e0 + v0], m_mesh.m_elems[e0 + v1]};
            bool add =true;
            for (int p=0; p<2; p++){
                if (std::find(itt->second.xxnodes.begin(), itt->second.xxnodes.end(), n[p]) ==
                    itt->second.xxnodes.end()) {
                    add=false;
                }
            }
            if (add){
                PEdge edge ( n[0], n[1], dom);
                getsurface->add_edge(edge, m_mesh.m_materials[matdom].domain());
            }
        }
        for (int ve=0; ve < 8; ++ve) {
            int n = m_mesh.m_elems[e0 + ve];
            bool add=true;
            if (std::find(itt->second.xxnodes.begin(), itt->second.xxnodes.end(), n) ==
                itt->second.xxnodes.end()) {
                add=false;
            }
            if (add) {
                PVertex gid ( n, dom);
                getsurface->add_vertex(gid, m_mesh.m_materials[matdom].domain());
            }
        }
        if (!mesh.m_surfaces.empty( )) {
            std::map<std::string,Surface*>::iterator st=mesh.m_surfaces.find(mesh.m_surf_matname[itmatname]);
            if (st!=mesh.m_surfaces.end()) {
                mesh.m_surfaces.at(mesh.m_surf_matname[itmatname])=getsurface;
            }
            else{
                mesh.m_surfaces.insert(std::pair<std::string,Surface*> (mesh.m_surf_matname[itmatname],getsurface));
            }
        }
        else {
            mesh.m_surfaces.insert(std::pair<std::string,Surface*> (mesh.m_surf_matname[itmatname],getsurface));
        }
    }
#endif
}


bool Mesh3DPart::is_border_element(int el)
{
    for(index_t k=m_mesh.m_xadj[el];k<m_mesh.m_xadj[el+1];++k) {
        index_t neighbour = m_mesh.m_adjncy[k];
        if (m_mesh.elem_part(neighbour)!=m_proc) {
            return true;
        }
    }
    return false;
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
        index_t face_id = it->second;
        if (!m_face_border[face_id]) continue;
        // We have a face shared between cpus
        const PFace& fc = it->first;
        std::set<index_t> elemset;
        m_mesh.get_neighbour_elements(4, fc.n, elemset);
        std::set<index_t>::const_iterator it;
        // We emit a neighbouring face for each element that don't belong
        // to this cpu (as it's a face there should be one and only one
        for(it=elemset.begin();it!=elemset.end();++it) {
            if (m_mesh.elem_part(*it)==m_proc) continue;
            handle_neighbouring_face(face_id, fc, *it);
        }
    }
}

void Mesh3DPart::compute_edge_communications()
{
    // This is similare to compute_face_communications above
    edge_map_t::const_iterator it;
    edge_map_t edge_to_id(m_edge_to_id);
    for(it=edge_to_id.begin();it!=edge_to_id.end();++it) {
        index_t edge_id = it->second;
        if (!m_edge_border[edge_id]) continue;
        // We have an edge shared between cpus
        const PEdge& ed = it->first;
        std::set<index_t> elemset;
        m_mesh.get_neighbour_elements(2, ed.n, elemset);
        std::set<index_t>::const_iterator it;
        for(it=elemset.begin();it!=elemset.end();++it) {
            if (m_mesh.elem_part(*it)==m_proc) continue;
            handle_neighbouring_edge(edge_id, ed, *it);
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
        index_t nvx = it->second;
        std::set<index_t> elemset;
        m_mesh.get_neighbour_elements(1, vx.n, elemset);
        std::set<index_t>::const_iterator it;
        for(it=elemset.begin();it!=elemset.end();++it) {
            if (m_mesh.elem_part(*it)==m_proc) continue;
            handle_neighbouring_vertex(nvx, vx, *it);
        }
    }
}
void Mesh3DPart::handle_neighbouring_face(index_t lnf, const PFace& fc, index_t el)
{
    int dom = m_mesh.get_elem_domain(el);
    index_t e0 = m_mesh.m_elems_offs[el];
    const index_t *eldef = &m_mesh.m_elems[e0];
    int source_proc = m_mesh.elem_part(el);
    for(int nf=0;nf<6;++nf) {
        index_t n[4];
        for(int i=0;i<4;++i) n[i] = eldef[RefFace[nf].v[i]];
        PFace nfc(n,dom);
        if (!nfc.eq_geom(fc)) continue;
        /// add local facet
        m_comm[source_proc].m_faces[fc]=lnf;
        if (dom!=fc.domain()) {
            // if neighbour's domain is different add a new facet
            index_t inf = add_facet(nfc, true);
            m_comm[source_proc].m_faces[nfc]=inf;
        }
    }
}

void Mesh3DPart::handle_neighbouring_edge(index_t lne, const PEdge& ed, index_t el)
{
    index_t e0 = m_mesh.m_elems_offs[el];
    const index_t *eldef = &m_mesh.m_elems[e0];
    int source_proc = m_mesh.elem_part(el);
    for(int ne=0;ne<12;++ne) {
        index_t n[2];
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
            index_t ne = add_edge(loced, true);
            m_comm[source_proc].m_edges[loced] = ne;
        }
    }
}

void Mesh3DPart::handle_neighbouring_vertex(index_t lnv, const PVertex& vx, index_t el)
{
    int source_proc = m_mesh.elem_part(el);
    int dommask = m_mesh.m_vertex_domains[vx.n[0]];
    for(int dom=0;dom<=DM_MAX;++dom) {
        if ((dommask & (1<<dom))==0) continue;
        PVertex locvx(vx.n[0], dom);
        index_t nv = add_vertex(locvx, true);
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
        if (itfound==m_face_to_id.end())
            continue;
        // A face from this surface belongs to this processor
        const PFace& fc = it->first;
        new_surf->add_face(fc, itfound->second);
        int dom = fc.domain();
        // Also creates the edges and vertices
        for(int k=0;k<4;++k) {
            index_t va = fc.n[k];
            index_t vb = fc.n[(k+1)%4];
            PEdge ed(va, vb, dom);
            index_t eid = get(m_edge_to_id, ed, invalid_index);
            assert(eid!=invalid_index);
            new_surf->add_edge(ed, eid);
            PVertex vx(va, dom);
            index_t vid = get(m_vertex_to_id, vx, invalid_index);
            assert(vid!=invalid_index);
            new_surf->add_vertex(vx, vid);
        }
    }
    // We keep all surfaces even if locally empty
    m_surfaces.push_back(new_surf);
}

index_t Mesh3DPart::add_facet(const PFace& facet, bool border)
{
    index_t nf;
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

index_t Mesh3DPart::add_edge(const PEdge& edge, bool border)
{
    index_t ne;
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

index_t Mesh3DPart::add_vertex(const PVertex& vertex, bool border)
{
    index_t nv;
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

index_t Mesh3DPart::add_node(index_t v0)
{
    index_t nv;
    node_id_map_t::iterator it;
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

void Mesh3DPart::shape8_local2global(double const vco[3][8],
                                     const double& xi, const double& eta, const double& zeta,
                                     double& x, double& y, double& z) const
{
    x = 0.125 * (vco[0][0]*(1-xi)*(1-eta)*(1-zeta) + vco[0][1]*(1+xi)*(1-eta)*(1-zeta) +
                 vco[0][2]*(1+xi)*(1+eta)*(1-zeta) + vco[0][3]*(1-xi)*(1+eta)*(1-zeta) +
                 vco[0][4]*(1-xi)*(1-eta)*(1+zeta) + vco[0][5]*(1+xi)*(1-eta)*(1+zeta) +
                 vco[0][6]*(1+xi)*(1+eta)*(1+zeta) + vco[0][7]*(1-xi)*(1+eta)*(1+zeta));
    y = 0.125 * (vco[1][0]*(1-xi)*(1-eta)*(1-zeta) + vco[1][1]*(1+xi)*(1-eta)*(1-zeta) +
                 vco[1][2]*(1+xi)*(1+eta)*(1-zeta) + vco[1][3]*(1-xi)*(1+eta)*(1-zeta) +
                 vco[1][4]*(1-xi)*(1-eta)*(1+zeta) + vco[1][5]*(1+xi)*(1-eta)*(1+zeta) +
                 vco[1][6]*(1+xi)*(1+eta)*(1+zeta) + vco[1][7]*(1-xi)*(1+eta)*(1+zeta));
    z = 0.125 * (vco[2][0]*(1-xi)*(1-eta)*(1-zeta) + vco[2][1]*(1+xi)*(1-eta)*(1-zeta) +
                 vco[2][2]*(1+xi)*(1+eta)*(1-zeta) + vco[2][3]*(1-xi)*(1+eta)*(1-zeta) +
                 vco[2][4]*(1-xi)*(1-eta)*(1+zeta) + vco[2][5]*(1+xi)*(1-eta)*(1+zeta) +
                 vco[2][6]*(1+xi)*(1+eta)*(1+zeta) + vco[2][7]*(1-xi)*(1+eta)*(1+zeta));
}

void Mesh3DPart::compute_gll()
{
    if (!m_cfg || !m_cfg->use_mirror) return;
    if (!m_cfg->mirror_impl_surf) return;

    calcul_gll(m_cfg->ngll, m_gll);
}

void Mesh3DPart::handle_mirror(index_t el)
{
    if (!m_cfg || !m_cfg->use_mirror) return;
    if (!m_cfg->mirror_impl_surf) return;

    double vco[3][8];
    index_t e0 = m_mesh.m_elems_offs[el];
    for(int vt=0;vt<8;++vt) {
        index_t gv = m_mesh.m_elems[e0 + vt];
        vco[0][vt] = m_mesh.m_xco[gv];
        vco[1][vt] = m_mesh.m_yco[gv];
        vco[2][vt] = m_mesh.m_zco[gv];
    }

    double f0 = 0.;
    bool init = false;
    bool sign_pos = false;
    bool sign_minus = false;
    std::vector<index_t> mirror_e;
    std::vector<index_t> mirror_ijk;
    std::vector<double>  mirror_xyz;
    for(int k=0;k<m_cfg->ngll;++k) {
        for(int j=0;j<m_cfg->ngll;++j) {
            for(int i=0;i<m_cfg->ngll;++i) {
                double x, y, z;
                shape8_local2global(vco, m_gll[i], m_gll[j], m_gll[k], x, y, z);

                double r = m_cfg->mirror_impl_surf_radius;
                double xc = m_cfg->mirror_impl_surf_center[0];
                double yc = m_cfg->mirror_impl_surf_center[1];
                double zc = m_cfg->mirror_impl_surf_center[2];

                double dx = x-xc;
                double dy = y-yc;
                double dz = z-zc;

                double f = r*r - dx*dx + dy*dy + dz*dz;
                if (f >= 0.) {
                    mirror_e.push_back(el);
                    mirror_ijk.push_back(i);
                    mirror_ijk.push_back(j);
                    mirror_ijk.push_back(k);
                    mirror_xyz.push_back(x);
                    mirror_xyz.push_back(y);
                    mirror_xyz.push_back(z);
                    sign_pos = true;
                } else {
                    sign_minus = true;
                }
            }
        }
    }

    if (sign_pos && sign_minus) {
        m_mirror_e.insert  (m_mirror_e.end(),   mirror_e.begin(),   mirror_e.end()  );
        m_mirror_ijk.insert(m_mirror_ijk.end(), mirror_ijk.begin(), mirror_ijk.end());
        m_mirror_ijk.insert(m_mirror_ijk.end(), mirror_ijk.begin(), mirror_ijk.end());
        m_mirror_ijk.insert(m_mirror_ijk.end(), mirror_ijk.begin(), mirror_ijk.end());
        m_mirror_xyz.insert(m_mirror_xyz.end(), mirror_xyz.begin(), mirror_xyz.end());
        m_mirror_xyz.insert(m_mirror_xyz.end(), mirror_xyz.begin(), mirror_xyz.end());
        m_mirror_xyz.insert(m_mirror_xyz.end(), mirror_xyz.begin(), mirror_xyz.end());
    }
}

void Mesh3DPart::handle_local_element(index_t el, bool is_border)
{
    index_t e0 = m_mesh.m_elems_offs[el];
    int dom = m_mesh.get_elem_domain(el);
    bool dom0 = false;
    surf_info_map_t::const_iterator it;
    it = m_mesh.surfelem.find(el);

    if ((dom0==false)&&(it==m_mesh.surfelem.end())){

        m_elems.push_back(el);
        // Assign all 6 faces
        for(int fc=0;fc<6;++fc) {
            index_t n[4];
            for(int p=0;p<4;++p) {
                n[p] = m_mesh.m_elems[e0 + RefFace[fc].v[p]];
            }
            PFace facet(n, dom);
            index_t nf = add_facet(facet, is_border);
            m_elems_faces.push_back(nf);
        }
        for(int ed=0;ed<12;++ed) {
            int v0 = RefEdge[ed][0];
            int v1 = RefEdge[ed][1];
            PEdge edge(m_mesh.m_elems[e0 + v0], m_mesh.m_elems[e0 + v1], dom);
            index_t ne = add_edge(edge, is_border);
            m_elems_edges.push_back(ne);
        }
        for(int vx=0;vx<8;++vx) {
            index_t gid = m_mesh.m_elems[e0 + vx];
            PVertex vertex(gid, dom);
            index_t vid = add_vertex(vertex, is_border);
            add_node(gid);
            m_elems_vertices.push_back(vid);
        }
        for(int vx=8;vx<m_mesh.nodes_per_elem();++vx) {
            index_t gid = m_mesh.m_elems[e0 + vx];
            add_node(gid);
        }
        handle_mirror(el);
    }
}


void Mesh3DPart::get_local_nodes(std::vector<double>& nodes) const
{
    node_id_map_t::const_iterator it;
    nodes.resize(3*m_nodes_to_id.size());
    for(it=m_nodes_to_id.begin();it!=m_nodes_to_id.end();++it) {
        int gid = it->first;
        int lid = it->second;
        nodes[3*lid+0] = m_mesh.m_xco[gid];
        nodes[3*lid+1] = m_mesh.m_yco[gid];
        nodes[3*lid+2] = m_mesh.m_zco[gid];
    }
}


void Mesh3DPart::get_local_elements(std::vector<loc_index_t>& elems) const
{
    int nctl_nodes = m_mesh.nodes_per_elem();
    elems.resize(m_elems.size()*nctl_nodes);

    for(size_t k=0;k<m_elems.size();++k) {
        int el = m_elems[k];
        int e0 = m_mesh.m_elems_offs[el];
        for(int n=0;n<nctl_nodes;++n) {
            index_t gid = m_mesh.m_elems[e0+n];
            index_t lid = get(m_nodes_to_id, gid, invalid_index);
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

void Mesh3DPart::get_local_mirrors(std::vector<int>& pmrrs) const
{
    pmrrs.resize(m_elems.size());
    size_t p=0;
    for(size_t k=0;k<m_elems.size();++k) {
        int el = m_elems[k];
        pmrrs[k] = m_mesh.m_mrrs[el];
    }
}

void Mesh3DPart::get_local_faces(std::vector<loc_index_t>& faces, std::vector<int>& doms) const
{
    faces.resize(m_face_to_id.size()*4);
    doms.resize(m_face_to_id.size());
    face_map_t::const_iterator it;
    index_t fc=0;
    for(it=m_face_to_id.begin();it!=m_face_to_id.end();++it) {
        fc = it->second;
        for(int p=0;p<4;++p) {
            index_t gid = it->first.n[p];
            index_t lid = get(m_nodes_to_id, gid, invalid_index);
            faces[4*fc + p] = lid;
            doms[fc] = it->first.domain();
        }
    }
}

void Mesh3DPart::get_local_edges(std::vector<loc_index_t>& edges, std::vector<int>& doms) const
{
    edges.resize(2*m_edge_to_id.size());
    doms.resize(m_edge_to_id.size());
    edge_map_t::const_iterator it;
    int ed=0;
    for(it=m_edge_to_id.begin();it!=m_edge_to_id.end();++it) {
        ed = it->second;
        for(int p=0;p<2;++p) {
            index_t gid = it->first.n[p];
            index_t lid = get(m_nodes_to_id, gid, invalid_index);
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

void Mesh3DPart::global_to_local_ids(std::vector<index_t>& ids) const
{
    node_id_map_t::const_iterator it, itend;
    itend = m_nodes_to_id.end();
    for(size_t k=0;k<ids.size();++k) {
        it = m_nodes_to_id.find(ids[k]);
        if (it!=itend)
            ids[k] = it->second;
    }
}

void Mesh3DPart::get_face_coupling(int d0, int d1, std::vector<index_t>& cpl, std::vector<int>& orient) const
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
        const PFace& fc0 = it0->first;
        const PFace& fc1 = it->first;
        index_t face_id0 = it0->second;
        index_t face_id1 = it->second;
        if (!fc1.eq_geom(fc0)) continue;
        // We have two faces equal with different domain,
        if (fc0.domain()!=d0) continue;
        if (fc1.domain()!=d1) continue;
        if (swapped) {
            cpl.push_back(face_id1);
            cpl.push_back(face_id0);
            orient.push_back(fc1.orient);
            orient.push_back(fc0.orient);
        } else {
            cpl.push_back(face_id0);
            cpl.push_back(face_id1);
            orient.push_back(fc0.orient);
            orient.push_back(fc1.orient);
        }
    }
}

void Mesh3DPart::get_edge_coupling(int d0, int d1, std::vector<index_t>& cpl) const
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

void Mesh3DPart::get_vertex_coupling(int d0, int d1, std::vector<index_t>& cpl) const
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
    h5h_create_attr(fid, "n_processors", (int)m_mesh.n_parts());
    h5h_create_attr(fid, "n_materials", (int)m_mesh.n_materials());
    h5h_create_attr(fid, "n_elements", (int)n_elems());
    h5h_create_attr(fid, "n_faces", (int)n_faces());
    h5h_create_attr(fid, "n_edges", (int)n_edges());
    h5h_create_attr(fid, "n_vertices", (int)n_vertices());
    h5h_create_attr(fid, "solid_fluid", false);
    h5h_create_attr(fid, "solid_fluid_loc", false);
    h5h_create_attr(fid, "all_fluid", false);
    h5h_create_attr(fid, "neumann_present", false);
    h5h_create_attr(fid, "n_surfaces", (int)m_mesh.n_surfaces("surface"));
    h5h_create_attr(fid, "curve", false);
}

static void convert_indexes(const std::vector<index_t>& src, std::vector<int>& dst)
{
    dst.resize(src.size());
    for(size_t k=0;k<src.size();++k) dst[k] = src[k];
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
    h5h_write_dset_2d(fid, "elements", n_elems(), m_mesh.nodes_per_elem(), tmpi.data());

    get_local_materials(tmpi,tmpi1);
    h5h_write_dset(fid, "material", n_elems(), tmpi.data());
    h5h_write_dset(fid, "domains", n_elems(), tmpi1.data());
    for(unsigned k=0;k<tmpi1.size();++k) {
        elem_doms[tmpi1[k]]++;
    }
    if (m_mesh.has_mrrs) {
        get_local_mirrors(tmpi);
        h5h_write_dset(fid, "mirror_pos", n_elems(), &tmpi[0]);
    }
    //
    convert_indexes(m_elems_faces, tmpi);
    h5h_write_dset_2d(fid, "faces", n_elems(), 6, tmpi.data());
    get_local_faces(tmpi, tmpi1);
    h5h_write_dset_2d(fid, "faces_def", n_faces(), 4, tmpi.data());
    h5h_write_dset(fid, "faces_dom", n_faces(), tmpi1.data());
    for(unsigned k=0;k<tmpi1.size();++k) {
        face_doms[tmpi1[k]]++;
    }
    //
    convert_indexes(m_elems_edges, tmpi);
    h5h_write_dset_2d(fid, "edges", n_elems(), 12, tmpi.data());
    get_local_edges(tmpi, tmpi1);
    h5h_write_dset_2d(fid, "edges_def", n_edges(), 2, tmpi.data());
    h5h_write_dset(fid, "edges_dom", n_edges(), tmpi1.data());
    for(unsigned k=0;k<tmpi1.size();++k) {
        edge_doms[tmpi1[k]]++;
    }
    convert_indexes(m_elems_vertices, tmpi);
    h5h_write_dset_2d(fid, "vertices", n_elems(), 8, tmpi.data());
    get_local_vertices_dom(tmpi);
    h5h_write_dset(fid, "vertices_dom", n_vertices(), tmpi.data());
    for(unsigned k=0;k<tmpi.size();++k) {
        vert_doms[tmpi[k]]++;
    }
    //
    printf("%04d : number of elements = (tot=%ld/fpml=%d/spml=%d/fl=%d/sol=%d)\n", m_proc, n_elems(),
           elem_doms[1],elem_doms[2],elem_doms[3],elem_doms[4]);
    printf("%04d : number of faces    = (tot=%ld/fpml=%d/spml=%d/fl=%d/sol=%d)\n", m_proc, n_faces(),
           face_doms[1],face_doms[2],face_doms[3],face_doms[4]);
    printf("%04d : number of edges    = (tot=%ld/fpml=%d/spml=%d/fl=%d/sol=%d)\n", m_proc, n_edges(),
           edge_doms[1],edge_doms[2],edge_doms[3],edge_doms[4]);
    printf("%04d : number of vertices = (tot=%ld/fpml=%d/spml=%d/fl=%d/sol=%d)\n", m_proc, n_vertices(),
           vert_doms[1],vert_doms[2],vert_doms[3],vert_doms[4]);
}

void Mesh3DPart::write_surface_dom(hid_t gid, const Surface* surf, const char* pfx, int dom)
{
    vector<int> data, orient, matdom, mat;
    //surf->surfelem_t =m_mesh.surfelem;

    char sface_data[100];
    char sface_orient[100];
//    char surf_matda[100];
    char sface_num [100];
    char sedge_data[100];
    char sedge_num [100];
    char svert_data[100];
    char svert_num [100];
    char svert_dom [100];
    char sedge_dom [100];
    char sface_dom [100];

    snprintf(sface_data, 100, "%s_faces", pfx);
    snprintf(sface_orient, 100, "%s_orient", pfx);
    snprintf(sface_num , 100, "n_%s_faces", pfx);
    snprintf(sedge_data, 100, "%s_edges", pfx);
    snprintf(sedge_num , 100, "n_%s_edges", pfx);
   // snprintf(surf_matda, 100, "%s_surf_mat", pfx);
    snprintf(svert_data, 100, "%s_vertices", pfx);
    snprintf(svert_num , 100, "n_%s_vertices", pfx);
    snprintf(svert_dom, 100, "%s_vertices_dom", pfx);
    snprintf(sedge_dom , 100, "%s_edges_dom", pfx);
    snprintf(sface_dom , 100, "%s_faces_dom", pfx);

    surf->get_faces_data(dom, data, orient, matdom);
    h5h_create_attr(gid, sface_num, (int)data.size());
    h5h_write_dset(gid, sface_data, data);
    h5h_write_dset(gid, sface_orient, orient);
    h5h_write_dset(gid, sface_dom, matdom);
    //h5h_write_dset(gid, surf_matda, mat);


    surf->get_edges_data(dom, data, orient, matdom);
    h5h_create_attr(gid, sedge_num, (int)data.size());
    h5h_write_dset(gid, sedge_data, data);
    h5h_write_dset(gid, sedge_dom, matdom);

    surf->get_vertices_data(dom, data, matdom);
    h5h_write_dset(gid, svert_data, data);
    h5h_create_attr(gid, svert_num, (int)data.size());
    h5h_write_dset(gid, svert_dom, matdom);
}


void Mesh3DPart::output_surface(hid_t fid, const Surface* surf)
{
    std::vector<int> faces, orient, asso_material;
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

    printf("%04d : number of elements = %ld\n", m_proc, n_elems());
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

    output_mirror();
}

void Mesh3DPart::output_mirror() const
{
    printf("%04d : number of mirror points = %ld\n", m_proc, m_mirror_xyz.size()/3);

    char fname[2048];
    snprintf(fname, sizeof(fname), "mirror.%04d.h5", m_proc);
    hid_t fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hid_t rid = H5Gcreate(fid, "Mirror", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5h_write_dset(rid, "E", m_mirror_e);
    //h5h_write_dset_2d(rid, "IJK", 3, m_mirror_ijk);
    //h5h_write_dset_2d(rid, "XYZ", 3, m_mirror_xyz);
    H5Gclose(rid);
    H5Fclose(fid);
}

void Mesh3DPart::output_comm(hid_t gid, const MeshPartComm& comm, int dest)
{
    std::vector<index_t> tdefs, tids;
    std::vector<int> tmpi, tdoms;
    h5h_create_attr(gid, "proc_dest", dest);
    h5h_create_attr(gid, "n_faces", (int)comm.n_faces());
    h5h_create_attr(gid, "n_edges", (int)comm.n_edges());
    h5h_create_attr(gid, "n_vertices", (int)comm.n_vertices());
    comm.get_faces(tdefs, tids, tdoms);
    convert_indexes(tids, tmpi);
    h5h_write_dset(gid, "faces", tmpi);
    if (m_mesh.debug) {
        global_to_local_ids(tdefs);
        convert_indexes(tdefs, tmpi);
        h5h_write_dset_2d(gid, "faces_def", 4, tmpi); // DEBUG
        h5h_write_dset(gid, "faces_dom", tdoms);
    }
    comm.get_edges(tdefs, tids, tdoms);
    convert_indexes(tids, tmpi);
    h5h_write_dset(gid, "edges", tmpi);
    if (m_mesh.debug) {
        global_to_local_ids(tdefs);
        convert_indexes(tdefs, tmpi);
        h5h_write_dset_2d(gid, "edges_def", 2, tmpi); // DEBUG
        h5h_write_dset(gid, "edges_dom", tdoms);
    }
    comm.get_vertices(tdefs, tids, tdoms);
    convert_indexes(tids, tmpi);
    h5h_write_dset(gid, "vertices", tmpi);
    if (m_mesh.debug) {
        global_to_local_ids(tdefs);
        convert_indexes(tdefs, tmpi);
        h5h_write_dset(gid, "vertices_defs", tmpi);
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
    vector<index_t> cpl;
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

    get_face_coupling(d0, d1, cpl, tmpo);
    convert_indexes(cpl, tmpi);
    h5h_create_attr(fid, sface_num, int(tmpi.size()/2) );
    h5h_write_dset_2d(fid, sface_data, tmpi.size()/2, 2, &tmpi[0]);
    h5h_write_dset_2d(fid, sface_orient, tmpo.size()/2, 2, &tmpo[0]);
    //
    get_edge_coupling(d0, d1, cpl);
    convert_indexes(cpl, tmpi);
    h5h_create_attr(fid, sedge_num, int(tmpi.size()/2) );
    h5h_write_dset_2d(fid, sedge_data, tmpi.size()/2, 2, &tmpi[0]);
    //
    get_vertex_coupling(d0, d1, cpl);
    convert_indexes(cpl, tmpi);
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
    output_xmf_mirror();
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
    fprintf(f, "      <Topology Type=\"Hexahedron\" NumberOfElements=\"%ld\">\n", n_elems());
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%ld 8\">\n", n_elems());
    fprintf(f, "mesh4spec.%04d.h5:/elements\n",m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Topology>\n");
    fprintf(f, "      <Geometry Type=\"XYZ\">\n");
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">\n", n_nodes());
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
    fprintf(f, "      <Topology Type=\"Quadrilateral\" NumberOfElements=\"%ld\">\n", n_faces());
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%ld 4\">\n", n_faces());
    fprintf(f, "mesh4spec.%04d.h5:/faces_def\n",m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Topology>\n");
    fprintf(f, "      <Geometry Type=\"XYZ\">\n");
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">\n", n_nodes());
    fprintf(f, "mesh4spec.%04d.h5:/local_nodes\n", m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Geometry>\n");
    output_int_scalar(f, 6, "FDom", "Cell", n_faces(), "/faces_dom");
    fprintf(f, "    </Grid>\n");
    output_xmf_footer(f);
    fclose(f);
}


void Mesh3DPart::output_xmf_mirror()
{
    char fname[2048];
    FILE* f;

    int n_faces_mirror = 0;
    for(unsigned k=0;k<m_surfaces.size();++k) {
        if (m_surfaces[k]->name()=="mirror") {
            n_faces_mirror = m_surfaces[k]->m_faces.size();
            break;
        }
    }

    snprintf(fname, sizeof(fname), "mesh4spec.%04d.mirror.xmf", m_proc);
    f = fopen(fname,"w");
    output_xmf_header(f);
    fprintf(f, "    <Grid Name=\"faces.%04d\">\n", m_proc);
    fprintf(f, "      <Topology Type=\"Quadrilateral\" NumberOfElements=\"%ld\">\n", n_faces());
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%ld 4\">\n", n_faces());
    fprintf(f, "mesh4spec.%04d.h5:/faces_def\n",m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Topology>\n");
    fprintf(f, "      <Geometry Type=\"XYZ\">\n");
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">\n", n_nodes());
    fprintf(f, "mesh4spec.%04d.h5:/local_nodes\n", m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Geometry>\n");
    output_int_scalar(f, 6, "FDom", "Cell", n_faces(), "/faces_dom");
    fprintf(f, "    </Grid>\n");
    fprintf(f, "    <Grid Name=\"surface.%04d\" GridType=\"Subset\" Section=\"DataItem\">\n", m_proc);
    fprintf(f, "      <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%d\">\n", n_faces_mirror);
    fprintf(f, "mesh4spec.%04d.h5:/Surfaces/mirror/sl_faces\n",m_proc);
    fprintf(f, "      </DataItem>\n");
    fprintf(f, "      <Grid Name=\"Target\" Reference=\"XML\">\n");
    fprintf(f, "/Xdmf/Domain/Grid[@Name=\"faces.%04d\"]\n", m_proc);
    fprintf(f, "      </Grid>\n");
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
    fprintf(f, "      <Topology Type=\"Polyline\" NumberOfElements=\"%ld\">\n", n_edges());
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Int\" Dimensions=\"%ld 2\">\n", n_edges());
    fprintf(f, "mesh4spec.%04d.h5:/edges_def\n",m_proc);
    fprintf(f, "        </DataItem>\n");
    fprintf(f, "      </Topology>\n");
    fprintf(f, "      <Geometry Type=\"XYZ\">\n");
    fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">\n", n_nodes());
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
        fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">\n", n_nodes());
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
        fprintf(f, "        <DataItem Format=\"HDF\" Datatype=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">\n", n_nodes());
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
