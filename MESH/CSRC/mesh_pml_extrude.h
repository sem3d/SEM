/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

// mesh_pml_extrude.h : add PML layers to an existing (imported) mesh by
// extruding boundary faces outward. Driven by an optional "pml.input" file.
#ifndef _MESH_PML_EXTRUDE_H_
#define _MESH_PML_EXTRUDE_H_

#include <string>

class Mesh3D;

// Directions, in the fixed processing order used by the extruder.
// x-/x+ then y-/y+ then z-/z+ so that corners/edges accumulate the flags of
// the sides already extruded (a W+S corner PML is created when the S pass
// meets the y-min face of a W column).
enum PmlSide { PML_XM=0, PML_XP, PML_YM, PML_YP, PML_ZM, PML_ZP, PML_NSIDES };

struct PmlSpec {
    int    n[PML_NSIDES];     // number of element layers per side (0 = off)
    double step[PML_NSIDES];  // extrusion step per side (<=0 -> auto from boundary element size)
    PmlSpec() { for(int k=0;k<PML_NSIDES;++k) { n[k]=0; step[k]=0.; } }
    bool any() const { for(int k=0;k<PML_NSIDES;++k) if(n[k]>0) return true; return false; }
};

// Parse pml.input. Returns true if the file exists and was read, false if it
// does not exist (caller then skips extrusion). Exits on malformed content.
bool read_pml_input(const std::string& fname, PmlSpec& spec);

// Extrude PML element layers onto mesh according to spec. Appends the derived
// PML materials to mesh.m_materials (so a subsequent write_materials emits them).
void extrude_pml(Mesh3D& mesh, const PmlSpec& spec);

#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
