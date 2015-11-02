// UNV specifications : http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse

#include <read_unv.hpp>
#include <iostream>
#include <fstream>
#include <sstream>

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

  vector<unsigned long long int> nodeID;
  for ( unsigned long long int i = 0; i < nbnodes; i++ )
  {
    unsigned long long int nID = 0; file >> nID; nID--; // nID: 1-based to 0-based
    if ( !file || nID >= nodes.size () ) { cerr << "read_unv_mesh_elems : bad element (" << elems.size () - 1 << ")" << endl; return 1; }
    nodeID.push_back ( nID );
  }
  int dim = 1; if ( type > 40 ) dim = 2; if ( type > 100 ) dim = 3; // Types are ordered : 1D < 40, 40 <= 2D < 100, 100 <= 3D
  elems.push_back ( elem { type, dim, nodeID, group { "", -1 } } );
  return 0;
}

int read_unv_mesh_groups ( ifstream & file, lsnodes & nodes, lselems & elems, lsgroups & gps )
{
  int dummy = 0; file >> dummy; if ( !file ) return 1; if ( dummy == -1 ) return -1; // Stop if block's -1 is reached

  unsigned long long int nbelems = 0; string gpname;
  file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> nbelems >> gpname;
  if ( !file ) { cerr << "read_unv_mesh_groups : bad group header (" << gpname << ")" << endl; return 1; }

  group gp { gpname, gps.size () }; bool added = false;
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
  if ( added ) gps.push_back ( gp );
  return 0;
}

int read_unv_mesh_block ( string const & fname, int const typeID, lsnodes & nodes, lselems & elems, lsgroups & groups )
{
  ifstream file ( fname.c_str () );
  int minus1 = 0; file >> minus1;
  while ( file )
  {
    int type = 0; file >> type;
    if ( minus1 == -1 && type > 0 ) // Read file block by block
    {
      int end = 0;
      do
      {
        if      ( type == typeID && typeID == 2411 ) { while ( end == 0 ) { end = read_unv_mesh_nodes  ( file, nodes );                } }
        else if ( type == typeID && typeID == 2412 ) { while ( end == 0 ) { end = read_unv_mesh_elems  ( file, nodes, elems );         } }
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

int read_unv_mesh ( string const & fname, lsnodes & nodes, lselems & elems )
{
  lsgroups groups;
  if ( read_unv_mesh_block ( fname, 2411, nodes, elems, groups ) != 0 ) return 1; // Nodes
  cout << "read_unv_mesh : " << nodes.size () << " nodes " << endl;
  if ( read_unv_mesh_block ( fname, 2412, nodes, elems, groups ) != 0 ) return 1; // Elements
  cout << "read_unv_mesh : " << elems.size () << " elements" << endl;
  if ( read_unv_mesh_block ( fname, 2467, nodes, elems, groups ) != 0 ) return 1; // Groups
  for ( unsigned int i = 0; i < groups.size (); i++ ) cout << "read_unv_mesh : group " << get<0> ( groups[i] ) << " has been found" << endl;
  return 0;
}

//int main ( int argc, char ** argv ) { if ( argc == 2 ) { lsnodes nodes; lselems elems; read_unv_mesh ( argv[1], nodes, elems ); } return 0; } // Debug read_unv_mesh
