#include <cstdio>
#include "mesh.h"
#include "readers.h"
#include "mesh_common.h"

class UNV_reader
{
public:
    void read_ideas_file(Mesh3D& mesh, const char* fname, int default_mat);
protected:
    int m_default_mat;
    std::vector<int> nodeNb;
    std::vector<int> elem_UNV;
    std::vector<int> mat_UNV;
};


void UNV_reader::read_ideas_file(Mesh3D& mesh, const char* fname, int default_mat)
{
    char buffer[2048];
    FILE* f = fopen(fname, "r");
    int desc;
    char* test

    m_default_mat = default_mat;

   do
   {
     test = fgets(buffer, 2048, f);
     sscanf(buffer, "%d",  &desc );

     if(desc == 2411 || desc == 781)
     {
         read_ideas_point_coords(mesh, f);
     }
     else if (desc == 2412 || desc == 780 )
     {
         read_connectivity(mesh, f);
     }
     else if (desc == 2477)
     {
         read_physical_volume(mesh, f)
     }
   } while (test != NULL)


   fclose(f);
}

void read_ideas_point_coords(Mesh3D& mesh, FILE* f)
{
    char buffer[2048];
    double x, y, z;  
    int i, a;
    int baseBlock, nodeNum;
    char* test

    i = 0;
    baseBlock = 2000;
    nodeNb.reserve(baseBlock);

    do
    {
      test = fgets(buffer, 2048, f);
      sscanf(buffer, "%d", a)
      if(a < 1)
      {
          break;
      }
      else
      {
          sscanf(buffer, "%d", &nodeNum)
 
          if(nodeNb.size() <= nodeNum)
          {
              nodeNb.resize(nodeNum+1, -1);
             
          }
          fgets(buffer, 2048, f);
          sscanf(buffer, "%f,%f,%f",  &x, &y, &z);
          nodeNb[nodeNum] = mesh.add_node(x, y, z);
      } while (test != NULL)
    }
}

void read_connectivity(Mesh3D& mesh, FILE* f)
{
    char buffer[2048];
    int a, nElem, typeElem, foo1, foo2, foo3, nNodes
    HexElem elem

    do
    {
      fgets(buffer, 2048, f);
      sscanf(buffer, "%d", a)
      if(a < 1)
      {
          break;
      }
      else
      {
          sscanf(buffer, "%d, %d, %d, %d, %d, %d", &nElem, &typeElem, &foo1, &foo2, &foo3, &nNodes)

          fgets(buffer, 2048, f)
          sscanf(buffer, "%d, %d, %d, %d, %d, %d", elem.v[0], elem.v[1], elem.v[2], elem.v[3], elem.v[4], elem.v[5], elem.v[6], elem.v[7] )
          add_elem(m_mat, elem)          
          

      }
    }    
}

void read_physical_volume(Mesh3D& mesh, FILE* f)
{
    char buffer[2048];
    int a, nElem, typeElem, foo1, foo2, foo3, nNodes
    HexElem elem

    do
    {
      fgets(buffer, 2048, f);
      sscanf(buffer, "%d", a)
      if(a < 1)
      {
          break;
      }
      else
      {
          sscanf(buffer, "%d, %d, %d, %d, %d, %d", &nElem, &typeElem, &foo1, &foo2, &foo3, &nNodes)

          fgets(buffer, 2048, f)
          sscanf(buffer, "%d, %d, %d, %d, %d, %d", elem.v[0], elem.v[1], elem.v[2], elem.v[3], elem.v[4], elem.v[5], elem.v[6], elem.v[7] )
          add_elem(m_mat, elem)          
          

      }
    }    
}


