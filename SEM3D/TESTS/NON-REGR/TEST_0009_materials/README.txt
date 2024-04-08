
This is an example of how to use extended material definitions in SEM :

- use material.spec instead of material.input to define materials

- Properties can be specified in an HDF5 file. The script create_material.py shows
  how to create a simple file with properties.

- Properties can be defined per-domain or be reused by several domains.

- Properties are defined on a cartesian (or polar) grid and linearly interpolated
  xMinGlob, xMaxGlob are the physical boundaries of the grid.

  A property can be used as long as a material domain (from the mesh) lies within the boundaries defined by
   xMinGlob, xMaxGlob.


