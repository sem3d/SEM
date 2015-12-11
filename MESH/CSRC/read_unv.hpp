// UNV specifications : http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse

#ifndef read_unv_h
#define read_unv_h 1

#include <string>
#include <vector>
#include <tuple>

using namespace std;

typedef tuple<string, int>                                     group;    // Group : name, group ID (0-based)
typedef tuple<double, double, double, group>                   node;     // Nodes : X, Y, Z, group
typedef tuple<int, int, vector<unsigned long long int>, group> elem;     // Elements : type, dim, node ID (0-based), group
typedef vector<group>                                          lsgroups;
typedef vector<node>                                           lsnodes;
typedef vector<elem>                                           lselems;

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
int read_unv_mesh ( string const & fpath, lsnodes & nodes, lselems & elems,
                    vector<int> * filterelems = NULL, vector<string> * filtergroups = NULL,
                    bool vtkformat = true );

#endif
