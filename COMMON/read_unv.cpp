// UNV specifications : http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse

#include <read_unv.hpp>
#include <iostream>      // cout, cerr
#include <fstream>       // ifstream
#include <sstream>       // stringstream
#include <unordered_map> // unordered_map
#include <algorithm>     // find

using namespace std;

int read_unv_mesh_nodes ( ifstream & file, lsnodes & nodes )
{
  int dummy = 0; file >> dummy; if ( !file ) return 1; if ( dummy == -1 ) return -1; // Stop if block's -1 is reached

  double x = 0., y = 0., z = 0.;
  file >> dummy >> dummy >> dummy >> x >> y >> z;
  if ( !file ) { cerr << "read_unv_mesh_nodes : bad node (" << nodes.size () - 1 << ")" << endl; return 1; }

  nodes.push_back ( node { x, y, z, group { "", -1 } } );
  return 0;
}

int read_unv_mesh_elems ( ifstream & file, lsnodes & nodes, lselems & elems )
{
  int dummy = 0; file >> dummy; if ( !file ) return 1; if ( dummy == -1 ) return -1; // Stop if block's -1 is reached

  unsigned long long int nbnodes = 0; int type = 0;
  file >> type >> dummy >> dummy >> dummy >> nbnodes;
  if ( type < 40 ) file >> dummy >> dummy >> dummy; // 1D element only (rod, beam : type < 40) : additional parameters
  if ( !file ) { cerr << "read_unv_mesh_elems : bad element header" << endl; return 1; }

  vector<unsigned long long int> nodeIDs;
  for ( unsigned long long int i = 0; i < nbnodes; i++ )
  {
    unsigned long long int nodeID = 0; file >> nodeID; nodeID--; // nodeID: 1-based to 0-based
    if ( !file || nodeID >= nodes.size () ) { cerr << "read_unv_mesh_elems : bad element (" << elems.size () - 1 << ")" << endl; return 1; }
    nodeIDs.push_back ( nodeID );
  }
  int dim = 1;
  if ( type > 40  ) dim = 2; if ( type > 110 ) dim = 3; // Types are ordered : 1D < 40, 2D >= 40, 3D >= 110
  if ( type > 130 ) dim = -1;                           // Does no more make sense if > 130
  elems.push_back ( elem { type, dim, nodeIDs, group { "", -1 } } );
  return 0;
}

int read_unv_mesh_groups ( ifstream & file, lsnodes & nodes, lselems & elems, lsgroups & groups )
{
  int dummy = 0; file >> dummy; if ( !file ) return 1; if ( dummy == -1 ) return -1; // Stop if block's -1 is reached

  unsigned long long int nbelems = 0; string gpname;
  file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> nbelems >> gpname;
  if ( !file ) { cerr << "read_unv_mesh_groups : bad group header (" << gpname << ")" << endl; return 1; }

  group gp { gpname, groups.size () }; bool added = false;
  for ( unsigned long long int i = 0; i < nbelems; i++ )
  {
    int type = 0; unsigned long long int tag = 0; file >> type >> tag >> dummy >> dummy; tag--; // tag: 1-based to 0-based
    if ( type == 7 ) // Group of nodes
    {
      if ( !file || tag >= nodes.size () ) { cerr << "read_unv_mesh_groups : bad node group (" << gpname << ", " << tag + 1 << ")" << endl; return 1; }
      get<3> ( nodes[tag] ) = gp; added = true;
    }
    if ( type == 8 ) // Group of elements
    {
      if ( !file || tag >= elems.size () ) { cerr << "read_unv_mesh_groups : bad element group (" << gpname << ", " << tag + 1 << ")" << endl; return 1; }
      get<3> ( elems[tag] ) = gp; added = true;
    }
  }
  if ( added ) groups.push_back ( gp );
  return 0;
}

int read_unv_mesh_block ( string const & fpath, int const typeID, lsnodes & nodes, lselems & elems, lsgroups & groups )
{
  ifstream file ( fpath.c_str () );
  int minus1 = 0; file >> minus1;
  while ( file )
  {
    int type = 0; file >> type;
    if ( minus1 == -1 && type > 0 ) // Read file block by block
    {
      int end = 0;
      do
      {
        if      ( type == typeID && typeID == 2411 ) { while ( end == 0 ) { end = read_unv_mesh_nodes  ( file, nodes                ); } }
        else if ( type == typeID && typeID == 2412 ) { while ( end == 0 ) { end = read_unv_mesh_elems  ( file, nodes, elems         ); } }
        else if ( type == typeID && typeID == 2467 ) { while ( end == 0 ) { end = read_unv_mesh_groups ( file, nodes, elems, groups ); } }
        else // Keep on reading file until target block is found
        {
          string line; getline ( file, line );
          stringstream sline; sline << line; sline >> end;
        }
      }
      while ( file && end != -1 );
      if ( end != -1 ) return 1; // KO
      type = end;                // OK, type is reset to "the next -1" as we reached the end of a block
    }
    minus1 = type;
  }
  return 0;
}

int filter_unv ( lsnodes & nodes, lselems & elems, vector<int> * filterelems, vector<string> * filtergroups )
{
  lselems pkelems; // Packed elements
  for ( unsigned int i = 0; i < elems.size (); i++ )
  {
    int type = get<0> ( elems[i] ); group gp = get<3> ( elems[i] ); string gpname = get<0> ( gp );
    if ( filterelems  && find ( filterelems  -> begin (), filterelems  -> end (), type   ) == filterelems  -> end () ) continue;
    if ( filtergroups && find ( filtergroups -> begin (), filtergroups -> end (), gpname ) == filtergroups -> end () ) continue;
    pkelems.push_back ( elems[i] );
  }
  elems = pkelems; // Replace with packed elements

  unordered_map<unsigned long long int, unsigned long long int> pkmap; // Packed map : non packed / packed ID
  lsnodes pknodes; ; // Packed nodes
  for ( unsigned int i = 0; i < elems.size (); i++ )
  {
    vector<unsigned long long int> nodeIDs = get<2> ( elems[i] );
    for ( unsigned int j = 0; j < nodeIDs.size (); j++ )
    {
      unsigned long long int nodeID = nodeIDs[j];
      if ( pkmap.find ( nodeID ) != pkmap.end () ) nodeIDs[j] = pkmap[nodeID];
      else { pknodes.push_back ( nodes[nodeID] ); nodeIDs[j] = pknodes.size () - 1; pkmap[nodeID] = nodeIDs[j]; }
    }
    get<2> ( elems[i] ) = nodeIDs; // Replace with packed node IDs
  }
  nodes = pknodes; // Replace with packed nodes
  return 0;
}

int transform_vtk ( lselems & elems )
{
  for ( unsigned int i = 0; i < elems.size (); i++ )
  {
    int type = get<0> ( elems[i] );
    vector<unsigned long long int> nodeIDs, vtkNodeIDs;

    if      ( type == 44 ) return 0; // Quad4 : OK
    else if ( type == 45 )           // Quad8
    {
      nodeIDs = get<2> ( elems[i] );
      for ( unsigned int j = 0; j < nodeIDs.size (); j += 2 ) vtkNodeIDs.push_back ( nodeIDs[j] ); // Main         nodes
      for ( unsigned int j = 1; j < nodeIDs.size (); j += 2 ) vtkNodeIDs.push_back ( nodeIDs[j] ); // Intermediate nodes
    }
    else { cerr << "transform_vtk : element type " << type << " not yet implemented" << endl; return 1; }

    get<2> ( elems[i] ) = vtkNodeIDs; // Replace node IDs with vtk node IDs
  }
  return 0;
}

int read_unv_mesh ( string const & fpath, lsnodes & nodes, lselems & elems, vector<int> * filterelems, vector<string> * filtergroups, bool vtkformat )
{
  lsgroups groups;
  if ( read_unv_mesh_block ( fpath, 2411, nodes, elems, groups ) != 0 ) return 1; // Nodes
  if ( read_unv_mesh_block ( fpath, 2412, nodes, elems, groups ) != 0 ) return 1; // Elements
  if ( read_unv_mesh_block ( fpath, 2467, nodes, elems, groups ) != 0 ) return 1; // Groups

  if ( ( filterelems || filtergroups ) && filter_unv    ( nodes, elems, filterelems, filtergroups ) != 0 ) return 1;
  if (   vtkformat                     && transform_vtk (        elems                            ) != 0 ) return 1;

  cout << "read_unv_mesh : " << nodes.size () << " nodes " << endl;
  cout << "read_unv_mesh : " << elems.size () << " elements" << endl;
  for ( unsigned int i = 0; i < groups.size (); i++ ) cout << "read_unv_mesh : group " << get<0> ( groups[i] ) << " has been found" << endl;
  return 0;
}

//int main ( int argc, char ** argv ) { if ( argc == 2 ) { lsnodes nodes; lselems elems; read_unv_mesh ( argv[1], nodes, elems ); } return 0; } // Debug read_unv_mesh
