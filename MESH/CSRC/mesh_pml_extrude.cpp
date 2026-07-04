/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include "mesh_pml_extrude.h"
#include "mesh.h"
#include "mesh_common.h"
#include "meshbase.h"
#include "material.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <vector>
#include <utility>

using namespace std;

// ---- pml.input parsing --------------------------------------------------
// Format (comments start with '#'):
//   x- 5
//   x+ 5 250.       # optional TOTAL PML thickness on this side (here 5 layers over 250 m)
//   z- 4
// A side is off unless listed. Missing 3rd value -> per-layer size = boundary element size.
static int side_from_token(const char* tok)
{
    if (!strcmp(tok,"x-")) return PML_XM;
    if (!strcmp(tok,"x+")) return PML_XP;
    if (!strcmp(tok,"y-")) return PML_YM;
    if (!strcmp(tok,"y+")) return PML_YP;
    if (!strcmp(tok,"z-")) return PML_ZM;
    if (!strcmp(tok,"z+")) return PML_ZP;
    return -1;
}

bool read_pml_input(const string& fname, PmlSpec& spec)
{
    FILE* f = fopen(fname.c_str(), "r");
    if (!f) return false;
    char* buffer=NULL; size_t linesize=0;
    while (true) {
        getData_line(&buffer, &linesize, f);
        if (!buffer || buffer[0]==0) break;
        char tok[64]; int n=0; double step=0.;
        int c = sscanf(buffer, "%63s %d %lf", tok, &n, &step);
        if (c<2) continue;
        int s = side_from_token(tok);
        if (s<0) { printf("ERR pml.input: unknown side '%s'\n", tok); exit(1); }
        if (n<0)  { printf("ERR pml.input: negative count for '%s'\n", tok); exit(1); }
        spec.n[s]    = n;
        spec.step[s] = (c>=3) ? step : 0.;
    }
    if (buffer) free(buffer);
    printf("Read pml.input: x-=%d x+=%d y-=%d y+=%d z-=%d z+=%d\n",
           spec.n[PML_XM], spec.n[PML_XP], spec.n[PML_YM],
           spec.n[PML_YP], spec.n[PML_ZM], spec.n[PML_ZP]);
    return true;
}

// ---- extrusion ----------------------------------------------------------
namespace {

// Per-material PML bookkeeping, indexed like mesh.m_materials.
struct MatInfo {
    int  base;                       // underlying non-PML material
    bool W,E,S,N,U,D;                // PML flags already set
    double xpos,xwidth,ypos,ywidth,zpos,zwidth;
    MatInfo():base(-1),W(false),E(false),S(false),N(false),U(false),D(false),
              xpos(0),xwidth(0),ypos(0),ywidth(0),zpos(0),zwidth(0) {}
};

struct SideGeom { int axis; double sign; }; // outward direction
static const SideGeom SIDE_GEOM[PML_NSIDES] = {
    {0,-1.}, {0,+1.}, {1,-1.}, {1,+1.}, {2,-1.}, {2,+1.}
};
static const int SIDE_REFFACE[PML_NSIDES] = { 4, 2, 1, 3, 0, 5 };

class PmlExtruder {
public:
    PmlExtruder(Mesh3D& m):mesh(m) { seed_materials(); }
    void run(const PmlSpec& spec);

private:
    Mesh3D& mesh;
    vector<MatInfo> minfo;
    // (original boundary node, layer) -> extruded node id
    map<pair<index_t,int>, index_t> newnode;

    double coord(index_t node, int axis) const {
        if (axis==0) return mesh.m_xco[node];
        if (axis==1) return mesh.m_yco[node];
        return mesh.m_zco[node];
    }
    void set_coord(double xyz[3], index_t node) const {
        xyz[0]=mesh.m_xco[node]; xyz[1]=mesh.m_yco[node]; xyz[2]=mesh.m_zco[node];
    }
    void seed_materials();
    int  get_or_make_pml(int src_mat, int side, double pos, double width, int axis);
    index_t extruded_node(index_t orig, int layer, int axis, double delta);
    void extrude_side(int side, int n, double step_override);
    void emit_hex(const index_t nodes8[8], int mat);
};

void PmlExtruder::seed_materials()
{
    minfo.resize(mesh.m_materials.size());
    for(size_t k=0;k<mesh.m_materials.size();++k) {
        const Material& mat = mesh.m_materials[k];
        MatInfo& mi = minfo[k];
        if (mat.is_pml()) {
            mi.base = (mat.associated_material>=0) ? mat.associated_material : (int)k;
            mi.W = mat.xwidth<0; mi.E = mat.xwidth>0;
            mi.S = mat.ywidth<0; mi.N = mat.ywidth>0;
            mi.D = mat.zwidth<0; mi.U = mat.zwidth>0;
            mi.xpos=mat.xpos; mi.xwidth=mat.xwidth;
            mi.ypos=mat.ypos; mi.ywidth=mat.ywidth;
            mi.zpos=mat.zpos; mi.zwidth=mat.zwidth;
        } else {
            mi.base = (int)k;
        }
    }
}

int PmlExtruder::get_or_make_pml(int src_mat, int side, double pos, double width, int axis)
{
    MatInfo si = minfo[src_mat];              // copy: base and accumulated flags/borders
    int base = si.base;
    bool W=si.W,E=si.E,S=si.S,N=si.N,U=si.U,D=si.D;
    bool* flag[PML_NSIDES] = {&W,&E,&S,&N,&D,&U}; // order matches PmlSide (ZM=D, ZP=U)
    if (*flag[side]) return src_mat;          // already a PML for this side
    *flag[side] = true;

    double xpos=si.xpos,xwidth=si.xwidth,ypos=si.ypos,ywidth=si.ywidth,zpos=si.zpos,zwidth=si.zwidth;
    if (axis==0) { xpos=pos; xwidth=width; }
    else if (axis==1) { ypos=pos; ywidth=width; }
    else { zpos=pos; zwidth=width; }

    Material& bm = mesh.m_materials[base];
    int idx = bm.pml_idx(W,E,S,N,U,D);
    int cached = bm.m_pml_num[idx];
    if (cached>=0) return cached;

    Material nm(mesh.m_materials[base]);
    if (nm.m_type==DM_SOLID_DG || nm.m_type==DM_SOLID_CG) nm.m_type = DM_SOLID_CG_PML;
    if (nm.m_type==DM_FLUID_DG || nm.m_type==DM_FLUID_CG) nm.m_type = DM_FLUID_CG_PML;
    nm.cinitial_type = nm.material_char();
    nm.set_pml_borders(xpos,xwidth,ypos,ywidth,zpos,zwidth);
    nm.associated_material = base;

    int newidx = mesh.m_materials.size();
    bm.m_pml_num[idx] = newidx;               // bm still valid (push_back below)
    mesh.m_materials.push_back(nm);

    MatInfo ni;
    ni.base=base; ni.W=W;ni.E=E;ni.S=S;ni.N=N;ni.U=U;ni.D=D;
    ni.xpos=xpos;ni.xwidth=xwidth;ni.ypos=ypos;ni.ywidth=ywidth;ni.zpos=zpos;ni.zwidth=zwidth;
    minfo.push_back(ni);
    printf("  PML material %d (base %d) flags W%d E%d S%d N%d U%d D%d\n",
           newidx, base, W,E,S,N,U,D);
    return newidx;
}

index_t PmlExtruder::extruded_node(index_t orig, int layer, int axis, double delta)
{
    pair<index_t,int> key(orig, layer);
    map<pair<index_t,int>,index_t>::iterator it = newnode.find(key);
    if (it!=newnode.end()) return it->second;
    double xyz[3]; set_coord(xyz, orig);
    xyz[axis] += delta*layer;
    index_t nid = mesh.add_node(xyz[0], xyz[1], xyz[2]);
    newnode[key] = nid;
    return nid;
}

// Order the 8 axis-aligned corners into the canonical SEM hex numbering
// (v0-3 = low-z face CCW-from-top, v4-7 = high-z face) so that later passes
// can find boundary faces of these elements via RefFace. Assumes the hex is
// axis-aligned (true for extruded PML layers): every corner is a distinct
// combination of {min,max} on each axis.
void PmlExtruder::emit_hex(const index_t nodes8[8], int mat)
{
    double mn[3], mx[3];
    for(int a=0;a<3;++a) { mn[a]=coord(nodes8[0],a); mx[a]=mn[a]; }
    for(int i=1;i<8;++i) for(int a=0;a<3;++a) {
        double c=coord(nodes8[i],a); if(c<mn[a])mn[a]=c; if(c>mx[a])mx[a]=c;
    }
    double mid[3]; for(int a=0;a<3;++a) mid[a]=0.5*(mn[a]+mx[a]);
    // slot by corner signature: bit0=xH, bit1=yH, bit2=zH -> SEM slot
    static const int slot_of[8] = {0,1,3,2,4,5,7,6};
    int v[8]; for(int k=0;k<8;++k) v[k]=-1;
    for(int i=0;i<8;++i) {
        int bx = coord(nodes8[i],0)>mid[0];
        int by = coord(nodes8[i],1)>mid[1];
        int bz = coord(nodes8[i],2)>mid[2];
        int slot = slot_of[bx | (by<<1) | (bz<<2)];
        v[slot] = nodes8[i];
    }
    Elem el(8);
    for(int k=0;k<8;++k) {
        if (v[k]<0) { printf("ERR: degenerate extruded hex (non axis-aligned corner)\n"); exit(1); }
        el.v[k]=v[k];
    }
    mesh.add_elem(mat, el);
}

void PmlExtruder::extrude_side(int side, int n, double step_override)
{
    if (n<=0) return;
    if (mesh.nodes_per_elem()!=8) {
        printf("ERR: PML extrusion only supports 8-node hexahedra (got %d control nodes)\n",
               mesh.nodes_per_elem());
        exit(1);
    }
    int axis = SIDE_GEOM[side].axis;
    double sgn = SIDE_GEOM[side].sign;
    newnode.clear(); // extruded nodes are per-side (a (node,layer) key means a
                     // different point on each side); within-side dedup keeps corners conforming

    // Current global bbox along axis
    double lo=coord(0,axis), hi=lo, ext=0.;
    for(size_t k=0;k<mesh.n_vertices();++k) {
        double c=coord((index_t)k,axis);
        if (c<lo) lo=c; if (c>hi) hi=c;
    }
    ext = hi-lo; if (ext<=0.) { printf("ERR: degenerate mesh along axis %d\n", axis); exit(1); }
    double plane = (sgn<0.) ? lo : hi;
    double tol = 1e-6*ext;

    // Collect boundary faces on this side; determine step.
    // The optional pml.input value is the TOTAL PML thickness on this side -> per-layer
    // step = total/n. When omitted, the boundary element size is used (conforming PML).
    struct BFace { index_t f[4]; int mat; };
    vector<BFace> faces;
    double step = (step_override>0.) ? step_override/n : 0.;
    double tstep=0.;
    size_t n0 = mesh.n_elems();
    for(size_t e=0;e<n0;++e) {
        index_t nodes[8]; mesh.get_elem_nodes((index_t)e, nodes);
        const FaceDesc& rf = RefFace[SIDE_REFFACE[side]];
        bool on=true;
        for(int k=0;k<4;++k) if (fabs(coord(nodes[rf.v[k]],axis)-plane) > tol) { on=false; break; }
        if (!on) continue;
        // element thickness along axis (for auto step / uniformity check)
        double emin=coord(nodes[0],axis), emax=emin;
        for(int k=1;k<8;++k){ double c=coord(nodes[k],axis); if(c<emin)emin=c; if(c>emax)emax=c; }
        double th=emax-emin;
        if (tstep<=0.) tstep=th;
        else if (step_override<=0. && fabs(th-tstep) > 1e-3*tstep) {
            printf("ERR: boundary elements on side %d have non-uniform size (%g vs %g). "
                   "Give an explicit total PML thickness in pml.input.\n", side, th, tstep);
            exit(1);
        }
        BFace bf;
        for(int k=0;k<4;++k) bf.f[k]=nodes[rf.v[k]];
        bf.mat = mesh.m_mat[e];
        faces.push_back(bf);
    }
    if (faces.empty()) { printf("WARNING: no boundary faces found for side %d, skipped\n", side); return; }
    if (step<=0.) step=tstep;
    double delta = sgn*step;
    double width = sgn*step*n;
    printf("Side %d: %zu boundary faces, %d layers, step=%g\n", side, faces.size(), n, step);

    for(size_t i=0;i<faces.size();++i) {
        const BFace& bf = faces[i];
        int mat = get_or_make_pml(bf.mat, side, plane, width, axis);
        index_t nodes8[8];
        for(int l=1;l<=n;++l) {
            for(int k=0;k<4;++k) {
                nodes8[k]   = (l==1) ? bf.f[k] : extruded_node(bf.f[k], l-1, axis, delta);
                nodes8[k+4] = extruded_node(bf.f[k], l, axis, delta);
            }
            emit_hex(nodes8, mat);
        }
    }
}

} // namespace

void PmlExtruder::run(const PmlSpec& spec)
{
    for(int side=0; side<PML_NSIDES; ++side) {
        extrude_side(side, spec.n[side], spec.step[side]);
    }
}

void extrude_pml(Mesh3D& mesh, const PmlSpec& spec)
{
    if (!spec.any()) return;
    printf("Extruding PML layers onto imported mesh...\n");
    PmlExtruder ex(mesh);
    ex.run(spec);
}

void append_pml_material_spec(const string& fname, Mesh3D& mesh, size_t first_new)
{
    FILE* f = fopen(fname.c_str(), "r");   // only touch material.spec if it exists
    if (!f) return;
    fclose(f);
    f = fopen(fname.c_str(), "a");
    if (!f) { printf("WARNING: cannot append to %s\n", fname.c_str()); return; }
    fprintf(f, "\n# PML materials added by pml.input extrusion\n");
    int appended=0;
    for(size_t k=first_new; k<mesh.m_materials.size(); ++k) {
        const Material& m = mesh.m_materials[k];
        if (!m.is_pml()) continue;
        int base = (m.associated_material>=0) ? m.associated_material : (int)k;
        const char* dom = (m.domain()==DM_FLUID_CG_PML) ? "fluidpml" : "solidpml";
        fprintf(f, "material %zu {\n    copy = %d;\n    domain = %s;\n};\n", k, base, dom);
        appended++;
    }
    fclose(f);
    printf("Appended %d PML material block(s) to %s\n", appended, fname.c_str());
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
