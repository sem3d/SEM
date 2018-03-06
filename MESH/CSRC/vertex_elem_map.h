/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

// vertex_elem_map.h : manages a map from vertex to elem number
#ifndef _VERTEX_ELEM_MAP_H_
#define _VERTEX_ELEM_MAP_H_

#include <vector>
#include <map>
#include <set>
#include "meshbase.h"

/// Manages a map of vertex points to elements connected to them
class VertexElemMap
{
public:
    VertexElemMap() {}

    void init(size_t nelements) {
        m_elem_linked_list.reserve(nelements*8);
        m_elem_linked_list.push_back( std::pair<index_t,index_t>(0,0) );
    }
    void add_link(index_t vid, index_t eid) {
        index_t list_head = m_vertex_to_first_elem[vid];
        index_t new_head = m_elem_linked_list.size();
        if (list_head==0) {
            m_vertex_to_first_elem[vid] = new_head;
            m_elem_linked_list.push_back(std::pair<index_t,index_t>(eid, 0));
            return ;
        }
        // Since we will fill the structure element by element, there is no
        // need to check if an element has already been seen.
        m_elem_linked_list.push_back(std::pair<index_t,index_t>(eid, list_head));
        m_vertex_to_first_elem[vid] = new_head;
    }

    void vertex_to_elements(index_t vid, std::vector<index_t>& elements) const {
        std::map<index_t,index_t>::const_iterator it = m_vertex_to_first_elem.find(vid);
        if (it==m_vertex_to_first_elem.end()) return;
        index_t head = it->second;
        while(head!=0) {
            elements.push_back(m_elem_linked_list[head].first);
            head = m_elem_linked_list[head].second;
        }
    }
    void vertex_to_elements(index_t vid, std::set<index_t>& elements) const {
        std::map<index_t,index_t>::const_iterator it = m_vertex_to_first_elem.find(vid);
        if (it==m_vertex_to_first_elem.end()) return;
        index_t head = it->second;
        while(head!=0) {
            elements.insert(m_elem_linked_list[head].first);
            head = m_elem_linked_list[head].second;
        }
    }
protected:
    /// maps a vertex to an element in m_elem_linked_list
    /// which is the head of a linked list ended by -1
    /// we don't use element 0, as 
    std::map<index_t,index_t> m_vertex_to_first_elem;
    /// We use the pair as (fisrt:eid, second:next)
    std::vector<std::pair<index_t,index_t> > m_elem_linked_list;
};

#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */

