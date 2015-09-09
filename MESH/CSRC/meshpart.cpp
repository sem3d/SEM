
#include <cstdio>
#include "mesh.h"
#include "meshpart.h"




/// Reference numbering for faces
static FaceDesc RefFace[6] = {
    { {0,1,2,3}, {0, 1, 2, 3} },
    { {0,1,5,4}, {0, 4, 5, 6} },
    { {1,2,6,5}, {1, 7, 8, 4} },
    { {3,2,6,7}, {2, 7, 9,10} },
    { {0,3,7,4}, {3,10,11, 6} },
    { {4,5,6,7}, {5, 8, 9,11} },
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
    printf("Created %d facets\n", m_face_to_face_id.size());
}

int Mesh3DPart::add_facet(int n[4], int dom)
{
    int nf;
    face_map_t::iterator it;
    Face facet(n, dom);
    it = m_face_to_face_id.find(facet);
    if (it==m_face_to_face_id.end()) {
        // New face
        nf = m_face_to_face_id.size();
        m_face_to_face_id[facet] = nf;
    } else {
        nf = it->second;
    }
    return nf;
}

void Mesh3DPart::handle_local_element(int el)
{
    int e0 = m_mesh.m_elems_offs[el];
    int dom = m_mesh.get_elem_domain(el);
    // Assign all 6 faces
    for(int fc=0;fc<6;++fc) {
        int n[4];
        for(int p=0;p<4;++p) {
            n[p] = m_mesh.m_elems[e0 + RefFace[fc].v[p]];
        }
        int nf = add_facet(n, dom);
    }
}

void Mesh3DPart::handle_neighbour_element(int el)
{
    // use adjacency map to tell if a neigbouring element belong to this proc
    // this works only if we used ncommon=1 in MeshToDual
    bool found = false;
    for(int k=m_mesh.m_xadj[el];k<m_mesh.m_xadj[el+1];++k) {
        int neighbour = m_mesh.m_adjncy[k];
        if (m_mesh.elem_part(neighbour)==m_proc) {
            found = true;
            break;
        }
    }
    if (!found) return;
    // We have a neighbouring element, 
}

void Mesh3DPart::output_mesh_part()
{
}

void Mesh3DPart::output_mesh_part_xmf()
{
}





/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
