/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include "mesh_grid.h"
#include "mesh.h"

using namespace std;

void RectMesh::read_params_old(FILE* fparam)
{
    char* buffer=NULL;
    size_t n=0;
    int pml_top, pml_bottom;

    getline(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &xmin);
    getline(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &xmax);
    getline(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &xstep);

    getline(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &ymin);
    getline(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &ymax);
    getline(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &ystep);

    getline(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &zmax);

    getline(&buffer, &n, fparam);
    sscanf(buffer, "%d", &nlayers);
    // Sanity checks
    if (nlayers<0 || nlayers>200) {
        printf("Check your parameter file, we read %d for nlayers\n", nlayers);
        exit(1);
    }
    thickness = (double*)malloc(nlayers*sizeof(double));
    nsteps    = (int*)   malloc(nlayers*sizeof(int));
    zmin = zmax;
    for(int k=0;k<nlayers;++k) {
        getline(&buffer, &n, fparam);
        sscanf(buffer, "%lf %d", &thickness[k], &nsteps[k]);
        zmin = zmin-thickness[k];
    }
    getline(&buffer, &n, fparam);
    sscanf(buffer, "%d", &has_pml);
    if (has_pml!=0 && has_pml!=1) {
        printf("Check your parameter file : we read has_pml=%d instead of 0 or 1\n", has_pml);
        exit(1);
    }
    if (has_pml) {
        pmls.N = true;
        pmls.S = true;
        pmls.E = true;
        pmls.W = true;
    } else {
        pmls.N = false;
        pmls.S = false;
        pmls.E = false;
        pmls.W = false;
    }
    getline(&buffer, &n, fparam);
    sscanf(buffer, "%d %d", &pml_top, &pml_bottom);
    if (has_pml && pml_top) {
        pmls.U = true;
    } else {
        pmls.U = false;
    }
    if (has_pml && pml_bottom) {
        pmls.D = true;
    } else {
        pmls.D = false;
    }
    getline(&buffer, &n, fparam);
    sscanf(buffer, "%d", &ngll_pml);
    getline(&buffer, &n, fparam);
    sscanf(buffer, "%d", &elem_shape);

}

int RectMesh::pointidx(int i, int j, int k)
{
    int idx =  i+j*(nelemx+1) + k*(nelemx+1)*(nelemy+1);
    //    printf("(%d / %d) (%d / %d)  (%d / %d)  = %d\n", i, nelemx, j, nelemy, k, 0, idx);
    return idx;
}

int RectMesh::get_mat(Mesh3D& mesh, int layer, bool W, bool E, bool S, bool N, bool U, bool D)
{
    if (!has_pml)
        return layer;
    W = pmls.W && W;
    E = pmls.E && E;
    S = pmls.S && S;
    N = pmls.N && N;
    U = pmls.U && U;
    D = pmls.D && D;

    bool pml = W||E||S||N||U||D;
    if (!pml)
        return layer;

    Material& mat = mesh.m_materials[layer];
    int new_idx = mat.pml_idx(W,E,S,N,U,D);
    int pml_mat = mat.m_pml_num[new_idx];
    if (pml_mat>=0) return pml_mat;

    Material new_mat(mat);
    if (new_mat.m_type==DM_SOLID) new_mat.m_type = DM_SOLID_PML;
    if (new_mat.m_type==DM_FLUID) new_mat.m_type = DM_FLUID_PML;
    double xw=0., yw=0., zw=0.;
    double xp=0., yp=0., zp=0.;

    if (W) { xw = -npml*xstep; xp = xmin0; }
    if (E) { xw =  npml*xstep; xp = xmax0; }
    if (S) { yw = -npml*ystep; yp = ymin0; }
    if (N) { yw =  npml*ystep; yp = ymax0; }
    if (U) { zw =  npml*thickness[0]/nsteps[0]; zp=zmax0; }
    if (D) { zw = -npml*thickness[nlayers-1]/nsteps[nlayers-1]; zp=zmin0; }

    new_mat.set_pml_borders(xp, xw, yp, yw, zp, zw);
    pml_mat = mesh.m_materials.size();
    printf("mat=%d layer = %d, idx=%d pmlnum=%d\n", pml_mat, layer, new_idx, mat.m_pml_num[new_idx]);
    printf("xp=%lf xw=%lf ; yp=%lf yw=%lf ; zp=%lf zw=%lf\n", xp, xw, yp, yw, zp, zw);
    mat.m_pml_num[new_idx] = pml_mat;
    mesh.m_materials.push_back(new_mat);
    return pml_mat;
}

void RectMesh::apply_pml_borders(int npml_)
{
    npml = npml_;
    xmax0 = xmax;
    xmin0 = ymin;
    ymax0 = ymax;
    ymin0 = ymin;
    zmax0 = zmax;
    zmin0 = zmin;
    if (pmls.E) { xmax+=npml*xstep; }
    if (pmls.W) { xmin-=npml*xstep; }
    if (pmls.N) { ymax+=npml*ystep; }
    if (pmls.S) { ymin-=npml*ystep; }
    if (pmls.U) {
        double zstep = thickness[0]/nsteps[0];
        zmax += npml*zstep;
        thickness[0] += npml*zstep;
        nsteps[0] += npml;
    }
    if (pmls.D) {
        int ll = nlayers-1;
        double zstep = thickness[ll]/nsteps[ll];
        thickness[ll] += npml*zstep;
        nsteps[ll] += npml;
    }
}

void RectMesh::init_rectangular_mesh(Mesh3D& mesh)
{
    const int npml=2;
    assert(xmin<xmax);
    assert(ymin<ymax);

    apply_pml_borders(npml);

    nelemx = int( (xmax-xmin)/xstep );
    nelemy = int( (ymax-ymin)/ystep );
    nelemz = 0;
    for(int k=0;k<nlayers;++k) {
        nelemz += nsteps[k];
    }

    printf("Creating grid mesh %d x %d x %d\n", nelemx, nelemy, nelemz );
    // Coordinates
    double layerzmax = zmax;
    int k0 = 0;
    int nnodes = 0;
    for(int nl=0;nl<nlayers;++nl) {
        double zmin = layerzmax - thickness[nl];
        double zstep = (layerzmax-zmin)/nsteps[nl];
        for(int k=k0;k<=nsteps[nl];++k) {
            for(int j=0;j<=nelemy;++j) {
                for(int i=0;i<=nelemx;++i) {
                    mesh.add_node( xmin + i*xstep,
                                   ymin + j*ystep,
                                   layerzmax - k*zstep );
                    nnodes++;
                }
            }
        }
        k0 = 1;
        layerzmax = zmin;
    }
    // Elements
    HexElem elem;
    int k=0;
    Surface* dirich = mesh.get_surface("dirichlet");
    for(int nl=0;nl<nlayers;++nl) {
        for(int kl=0;kl<nsteps[nl];++kl) {
            for(int j=0;j<nelemy;++j) {
                for(int i=0;i<nelemx;++i) {
                    elem.v[0] = pointidx(i  ,j  ,k+1); // bottom layer
                    elem.v[1] = pointidx(i+1,j  ,k+1);
                    elem.v[2] = pointidx(i+1,j+1,k+1);
                    elem.v[3] = pointidx(i  ,j+1,k+1);
                    elem.v[4] = pointidx(i  ,j  ,k  ); // top layer
                    elem.v[5] = pointidx(i+1,j  ,k  );
                    elem.v[6] = pointidx(i+1,j+1,k  );
                    elem.v[7] = pointidx(i  ,j+1,k  );
                    for(int in=0;in<8;++in) {
                        assert(elem.v[in]<nnodes);
                        assert(elem.v[in]>=0);
                    }
                    bool W = (i<npml);
                    bool E = (i>(nelemx-npml-1));
                    bool S = (j<npml);
                    bool N = (j>(nelemy-npml-1));
                    bool U = (k<npml && nl==0);
                    bool D = (k>(nelemz-npml-1) && nl==(nlayers-1));
                    bool LW = (i<1);
                    bool LE = (i>(nelemx-2));
                    bool LS = (j<1);
                    bool LN = (j>(nelemy-2));
                    bool LU = (k<1 && nl==0);
                    bool LD = (k>(nelemz-2) && nl==(nlayers-1));
                    int mat = get_mat(mesh, nl, W,E,S,N,U,D);
                    mesh.add_elem(mat, elem);
                    // Check for free fluid surface and free pml surface
                    int dom = mesh.m_materials[mat].domain();
                    if (dom==DM_FLUID) {
                        if (U) emit_free_face(dirich, dom, elem, 5);
                    }
                    if (dom==DM_FLUID_PML) {
                        emit_free_face(dirich, dom, elem, W,E,S,N,U,D);
                    }
                    if (dom==DM_SOLID_PML) {
                        double x_dir = mesh.m_materials[mat].xwidth;
                        double y_dir = mesh.m_materials[mat].ywidth;
                        double z_dir = mesh.m_materials[mat].zwidth;

                        emit_free_face(dirich, dom, elem,
                                       LW&&(x_dir<0), LE&&(x_dir>0),
                                       LS&&(y_dir<0), LN&&(y_dir>0),
                                       LU&&(z_dir>0), LD&&(z_dir<0));
                    }
                }
            }
            k++;
        }
    }
}

void RectMesh::emit_free_face(Surface* surf, int dom, const HexElem& elem,
                              bool W, bool E, bool S, bool N, bool U, bool D)
{
    if (W) emit_free_face(surf, dom, elem, 4);
    if (E) emit_free_face(surf, dom, elem, 2);
    if (S) emit_free_face(surf, dom, elem, 1);
    if (N) emit_free_face(surf, dom, elem, 3);
    if (U) emit_free_face(surf, dom, elem, 5);
    if (D) emit_free_face(surf, dom, elem, 0);
}

void RectMesh::emit_free_face(Surface* surf, int dom, const HexElem& elem, int facenum)
{
    int n[4];
    for(int k=0;k<4;++k) {
        n[k] = elem.v[RefFace[facenum].v[k]];
    }
    PFace fc(n, dom);
    surf->add_face(fc, 0);
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
