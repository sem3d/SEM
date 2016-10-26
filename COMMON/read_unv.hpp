// UNV specifications : http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse

#ifndef READ_UNV_HPP
#define READ_UNV_HPP

#include <string>
#include <vector>
#include <tuple>
#include <cstdint>

typedef std::tuple<std::string, int>                        group; // Group : name, group ID (0-based)
typedef std::tuple<double, double, double, group>           node;  // Nodes : X, Y, Z, group
typedef std::tuple<int, int, std::vector<uint64_t>, group>  elem;  // Elements : type, dim, node ID (0-based), group
typedef std::vector<group>                                  lsgroups;
typedef std::vector<node>                                   lsnodes;
typedef std::vector<elem>                                   lselems;

/*
 * read_unv_mesh : read mesh from unv file.
 * fpath        = path                  to   the unv file.
 * lsnodes      = list of nodes    read from the unv file.
 * lselems      = list of elements read from the unv file.
 * filterelems  = filter designed to include only elements of specific types.
 * filtergroups = filter designed to include only groups   of specific names.
 * vtkformat    = set returned elements to vtk format.
 *
 * If filters are not used, all  the nodes / elements contained in the unv file are returned.
 * If filters are     used, only the nodes / elements matching the filter       are returned (= a subset of the whole file).
 *
 * Note :
 * Usually, unv exported files are messy (they contain a lot more than one needs to use).
 * Filters enable to shrink the data set extracted from the unv file to a more convenient subset.
 */
int read_unv_mesh ( std::string const & fpath, lsnodes & nodes, lselems & elems,
                    std::vector<int> * filterelems = NULL,
                    std::vector<std::string> * filtergroups = NULL,
                    bool vtkformat = true );

#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
