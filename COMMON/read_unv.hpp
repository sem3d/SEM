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

int read_unv_mesh ( string const & fname, lsnodes & nodes, lselems & elems );

#endif
