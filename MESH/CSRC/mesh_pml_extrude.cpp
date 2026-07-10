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
#include <cassert>
#include <map>
#include <vector>
#include <utility>

using namespace std;

// ---- pml.input parsing --------------------------------------------------
// Format (comments start with '#'):  <side> [thickness] [n_elements] [ratio]
//   x- 250.         # 250 m thick, 1 element
//   x+ 250. 5       # 250 m thick, split into 5 elements
//   y- 250. 5 1.3   # + geometric grading 1.3 (finer near the domain, coarser outward)
//   z-              # thickness omitted -> one boundary element (conforming), warns
// A side is off unless listed. Thickness is the primary knob; if given it must be > 0.
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
        char tok[64]; int n=1; double thick=0., ratio=1.;
        if (sscanf(buffer, "%63s", tok)==1 && !strcmp(tok,"pmlparams")) {
            sscanf(buffer, "%*s %d %lf", &spec.npow, &spec.Rc);
            continue;
        }
        // <side> [thickness] [n_elements] [ratio]
        //   thickness  : TOTAL PML thickness on this side (the primary knob).
        //                Omitted -> one boundary element (conforming), with a warning
        //                nudge, since thickness matters more than the element count.
        //   n_elements : layers to split the thickness into (default 1).
        //   ratio      : geometric grading (default 1 = uniform; >1 = finer near domain).
        int c = sscanf(buffer, "%63s %lf %d %lf", tok, &thick, &n, &ratio);
        if (c<1) continue;
        int s = side_from_token(tok);
        if (s<0) { printf("ERR pml.input: unknown side '%s'\n", tok); exit(1); }
        if (c>=2 && thick<=0.) { printf("ERR pml.input: PML thickness for '%s' must be > 0 (got %g)\n", tok, thick); exit(1); }
        if (c>=3 && n<1)       { printf("ERR pml.input: element count for '%s' must be >= 1 (got %d)\n", tok, n); exit(1); }
        if (c>=4 && ratio<=0.) { printf("ERR pml.input: grading ratio for '%s' must be > 0 (got %g)\n", tok, ratio); exit(1); }
        if (c<2) printf("WARNING pml.input: no thickness for '%s'; using one boundary element. "
                        "Set a thickness -- it matters more than the element count.\n", tok);
        spec.n[s]     = (c>=3) ? n : 1;
        spec.step[s]  = (c>=2) ? thick : 0.; // 0 -> auto (one boundary element) in extrude_side
        spec.ratio[s] = (c>=4) ? ratio : 1.;
    }
    if (buffer) free(buffer);
    printf("Read pml.input: x-=%d x+=%d y-=%d y+=%d z-=%d z+=%d (npow=%d Rc=%g)\n",
           spec.n[PML_XM], spec.n[PML_XP], spec.n[PML_YM],
           spec.n[PML_YP], spec.n[PML_ZM], spec.n[PML_ZP], spec.npow, spec.Rc);
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
    int npow;
    double apow;
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
    index_t extruded_node(index_t orig, int layer, int axis, double offset);
    void extrude_side(int side, int n, double step_override, double ratio);
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
    nm.m_type = pml_domain_for(mesh.m_materials[base]); // solid/fluid PML from any base
    nm.cinitial_type = nm.material_char();
    nm.set_pml_borders(xpos,xwidth,ypos,ywidth,zpos,zwidth);
    nm.associated_material = base;
    nm.npow = npow;
    nm.apow = apow;

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

index_t PmlExtruder::extruded_node(index_t orig, int layer, int axis, double offset)
{
    pair<index_t,int> key(orig, layer);
    map<pair<index_t,int>,index_t>::iterator it = newnode.find(key);
    if (it!=newnode.end()) return it->second;
    double xyz[3]; set_coord(xyz, orig);
    xyz[axis] += offset; // absolute (already signed) offset of this layer plane
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

void PmlExtruder::extrude_side(int side, int n, double step_override, double ratio)
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

    // Total PML thickness on this side: explicit if given, else n boundary elements.
    double total = (step_override>0.) ? step_override : tstep*n;
    // Per-layer cumulative offsets off[0..n] (off[0]=0, off[n]=total).
    // Uniform when ratio==1; otherwise a geometric progression h_L = h1*ratio^(L-1),
    // with h1 = total*(r-1)/(r^n-1) so the layers grade from fine (near the domain)
    // to coarse (outer) while summing exactly to total.
    vector<double> off(n+1, 0.);
    if (fabs(ratio-1.) < 1e-9) {
        for(int l=1;l<=n;++l) off[l] = total*l/n;
    } else {
        double rn = pow(ratio, n);
        double h1 = total*(ratio-1.)/(rn-1.);
        double h = h1;
        for(int l=1;l<=n;++l) { off[l] = off[l-1] + h; h *= ratio; }
    }
    assert(fabs(off[n]-total) < 1e-9*total); // geometric series must reconstruct the total
    double width = sgn*off[n];
    printf("Side %d: %zu boundary faces, %d layers, total=%g ratio=%g (h1=%g h_n=%g)\n",
           side, faces.size(), n, total, ratio, off[1], off[n]-off[n-1]);

    for(size_t i=0;i<faces.size();++i) {
        const BFace& bf = faces[i];
        int mat = get_or_make_pml(bf.mat, side, plane, width, axis);
        index_t nodes8[8];
        for(int l=1;l<=n;++l) {
            for(int k=0;k<4;++k) {
                nodes8[k]   = (l==1) ? bf.f[k] : extruded_node(bf.f[k], l-1, axis, sgn*off[l-1]);
                nodes8[k+4] = extruded_node(bf.f[k], l, axis, sgn*off[l]);
            }
            emit_hex(nodes8, mat);
        }
    }
}

} // namespace

void PmlExtruder::run(const PmlSpec& spec)
{
    npow = spec.npow;
    apow = pml_apow_from_rc(npow, spec.Rc);
    for(int side=0; side<PML_NSIDES; ++side) {
        extrude_side(side, spec.n[side], spec.step[side], spec.ratio[side]);
    }
}

void extrude_pml(Mesh3D& mesh, const PmlSpec& spec)
{
    if (!spec.any()) return;
    printf("Extruding PML layers onto imported mesh...\n");
    PmlExtruder ex(mesh);
    ex.run(spec);
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
