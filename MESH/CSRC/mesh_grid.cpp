/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "mesh_grid.h"
#include "mesh.h"
#include "mesh_common.h"
#include "read_input.h"

using namespace std;

void RectMesh::read_params_old(FILE* fparam, sem_config_t* config)
{
    char* buffer=NULL;
    size_t n=0;

    //getline(&buffer, &n, fparam);
    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &xmin);
    //getline(&buffer, &n, fparam);
    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &xmax);
    //getline(&buffer, &n, fparam);
    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &xstep);

    //getline(&buffer, &n, fparam);
    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &ymin);
    //getline(&buffer, &n, fparam);
    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &ymax);
    //getline(&buffer, &n, fparam);
    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &ystep);

    //getline(&buffer, &n, fparam);
    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%lf", &zmax);

    //getline(&buffer, &n, fparam);
    getData_line(&buffer, &n, fparam);
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
        //getline(&buffer, &n, fparam);
    	getData_line(&buffer, &n, fparam);
    	sscanf(buffer, "%lf %d", &thickness[k], &nsteps[k]);
        zmin = zmin-thickness[k];
    }
    //getline(&buffer, &n, fparam);
    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%d", &has_pml);
    if (has_pml<0) {
        printf("Check your parameter file : we read has_pml=%d instead of 0 or 1\n", has_pml);
        exit(1);
    }
    if (has_pml>1) {
        printf("Using %d layers of PML\n", has_pml);
    }
    npml = has_pml;
    int pml_top, pml_bottom, pml_N, pml_S, pml_E, pml_W, flag;

    pmls.N = false;
    pmls.S = false;
    pmls.E = false;
    pmls.W = false;
    pmls.U = false;
    pmls.D = false;

    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%d", &pml_top);
    if (pml_top==2) {
        sscanf(buffer, "%d %d %d %d %d %d %d", &flag,
               &pml_top, &pml_bottom, &pml_N, &pml_S, &pml_E, &pml_W);
    } else {
        sscanf(buffer, "%d %d", &pml_top, &pml_bottom);
        pml_N = pml_S = pml_E = pml_W = 1;
    }

    if (has_pml) {
        if (pml_top)    pmls.U = true;
        if (pml_bottom) pmls.D = true;
        if (pml_N)      pmls.N = true;
        if (pml_S)      pmls.S = true;
        if (pml_E)      pmls.E = true;
        if (pml_W)      pmls.W = true;
    }

    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%d", &ngll_pml);

    getData_line(&buffer, &n, fparam);
    sscanf(buffer, "%d", &elem_shape);
    switch(elem_shape) {
    case 1:
    case 8:
        elem_shape = 8;
        break;
    case 2:
    case 27:
        elem_shape = 27;
        break;
    default:
        printf("Please use 1 or 8 for linear elements, or 2 or 27 for quadratic elements\n");
        exit(1);
    }

    // spherical shell top transf
    getData_line(&buffer, &n, fparam);
    has_sph=0;
    sscanf(buffer, "%d", &has_sph);

    // mirrors options
    getData_line(&buffer, &n, fparam);
    has_mirrors=0;
    sscanf(buffer, "%d", &has_mirrors);

    if (has_mirrors<0 || has_mirrors>2) {
        printf("Check your parameter file : we read has_mirros=%d instead of 0 or 1\n", has_mirrors);
        exit(1);
    }

    int m_W, m_E, m_S, m_N, m_D, m_U;
    if (has_mirrors == 1) {
        getData_line(&buffer, &n, fparam);
        sscanf(buffer, "%d %d %d %d %d %d", &m_W, &m_E, &m_S, &m_N, &m_D, &m_U);
        mrrs.W = mrrs.E = mrrs.S = mrrs.N = mrrs.D = mrrs.U = true;
        mrrs.nW = m_W;
        mrrs.nE = m_E;
        mrrs.nS = m_S;
        mrrs.nN = m_N;
        mrrs.nD = m_D;
        mrrs.nU = m_U;
        if (pml_top<1) {
            mrrs.U = false;
            mrrs.nU = -1;
        }
    }

    if (has_mirrors == 2) { // Read input.spec
        if (config) {
            int err;
            read_sem_config(config, 0, 3, "input.spec", &err);
            dump_config(config);
            printf("\n");
        }
    }
}

int RectMesh::pointidx(int i, int j, int k)
{
    int idx =  i+j*(nelemx+1) + k*(nelemx+1)*(nelemy+1);
    return idx;
}

int RectMesh::pointidx27(int i, int j, int k)
{
    int idx =  i+j*(2*nelemx+1) + k*(2*nelemx+1)*(2*nelemy+1);
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
    new_mat.cinitial_type = new_mat.material_char();
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
    new_mat.associated_material = layer;
    printf("mat=%d layer = %d, idx=%d pmlnum=%d\n", pml_mat, layer, new_idx, mat.m_pml_num[new_idx]);
    printf("xp=%lf xw=%lf ; yp=%lf yw=%lf ; zp=%lf zw=%lf\n", xp, xw, yp, yw, zp, zw);
    mat.m_pml_num[new_idx] = pml_mat;
    mesh.m_materials.push_back(new_mat);
    return pml_mat;
}

void RectMesh::apply_pml_borders()
{
    xmax0 = xmax;
    xmin0 = xmin;
    ymax0 = ymax;
    ymin0 = ymin;
    zmax0 = zmax;
    zmin0 = zmin;
    if (pmls.E) { xmax+=npml*xstep;}
    if (pmls.W) { xmin-=npml*xstep;}
    if (pmls.N) { ymax+=npml*ystep;}
    if (pmls.S) { ymin-=npml*ystep;}
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

    if (has_mirrors) {
        if (pmls.E) {mrrs.nE += npml;}
        if (pmls.W) {mrrs.nW += npml;}
        if (pmls.N) {mrrs.nN += npml;}
        if (pmls.S) {mrrs.nS += npml;}
        if (pmls.U) {mrrs.nU += npml;}
        if (pmls.D) {mrrs.nD += npml;}
    }
}

int RectMesh::create_linear_grid_nodes(Mesh3D& mesh)
{
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
    return nnodes;
}

int RectMesh::create_quadratic_grid_nodes(Mesh3D& mesh)
{
    // Coordinates
    double layerzmax = zmax;

    int k0 = 0;
    int nnodes = 0;
    int npx = 2*nelemx+1;
    int npy = 2*nelemy+1;

    for(int nl=0;nl<nlayers;++nl) {
        double zmin = layerzmax - thickness[nl];
        double zstep = (layerzmax-zmin)/nsteps[nl];
        for(int k=k0;k<2*nsteps[nl]+1;++k) {
            for(int j=0;j<npy;++j) {
                for(int i=0;i<npx;++i) {
                    mesh.add_node( xmin + i*0.5*xstep,
                                   ymin + j*0.5*ystep,
                                   layerzmax - k*0.5*zstep );
                    nnodes++;
                }
            }
        }
        k0 = 1;
        layerzmax = zmin;
    }
    return nnodes;
}

int RectMesh::create_linear_grid_nodes_sph(Mesh3D& mesh)
{
    // Coordinates
    double xloc,yloc,zloc,wloc;
    double layerzmax = zmax;
    double zmin = layerzmax;
    int k0 = 0;
    int nnodes = 0;
    for(int nl=0;nl<nlayers;++nl) {
        zmin -= thickness[nl];
    }
    for(int nl=0;nl<nlayers;++nl) {
        double layerzmin = layerzmax - thickness[nl];
        double zstep = (layerzmax-layerzmin)/nsteps[nl];
        for(int k=k0;k<=nsteps[nl];++k) {
            for(int j=0;j<=nelemy;++j) {
                for(int i=0;i<=nelemx;++i) {
                    xloc = xmin+i*xstep;
                    yloc = ymin+j*ystep;
                    zloc = sqrt(zmax*zmax-xloc*xloc-yloc*yloc);
                    wloc = (layerzmax-k*zstep-zmin)/(zmax-zmin);
                    zloc = wloc*zloc+(1.-wloc)*zmin;
                    mesh.add_node(xloc,yloc,zloc);
                    nnodes++;
                }
            }
        }
        k0 = 1;
        layerzmax = layerzmin;
    }
    return nnodes;
}

void RectMesh::create_linear_element(Elem& elem, int i, int j, int k)
{
    elem.v[0] = pointidx(i  ,j  ,k+1); // bottom layer
    elem.v[1] = pointidx(i+1,j  ,k+1);
    elem.v[2] = pointidx(i+1,j+1,k+1);
    elem.v[3] = pointidx(i  ,j+1,k+1);
    elem.v[4] = pointidx(i  ,j  ,k  ); // top layer
    elem.v[5] = pointidx(i+1,j  ,k  );
    elem.v[6] = pointidx(i+1,j+1,k  );
    elem.v[7] = pointidx(i  ,j+1,k  );
}

void RectMesh::create_quadratic_element(Elem& elem, int i, int j, int k)
{
    elem.v[ 0] = pointidx27(2*i  , 2*j  , 2*k+2); // bottom layer
    elem.v[ 1] = pointidx27(2*i+2, 2*j  , 2*k+2);
    elem.v[ 2] = pointidx27(2*i+2, 2*j+2, 2*k+2);
    elem.v[ 3] = pointidx27(2*i  , 2*j+2, 2*k+2);
    elem.v[ 4] = pointidx27(2*i  , 2*j  , 2*k  ); // top layer
    elem.v[ 5] = pointidx27(2*i+2, 2*j  , 2*k  );
    elem.v[ 6] = pointidx27(2*i+2, 2*j+2, 2*k  );
    elem.v[ 7] = pointidx27(2*i  , 2*j+2, 2*k  );

    elem.v[ 8] = pointidx27(2*i+1, 2*j+0, 2*k+2); // mid segment bottom layer
    elem.v[ 9] = pointidx27(2*i+2, 2*j+1, 2*k+2);
    elem.v[10] = pointidx27(2*i+1, 2*j+2, 2*k+2);
    elem.v[11] = pointidx27(2*i+0, 2*j+1, 2*k+2);

    elem.v[12] = pointidx27(2*i  , 2*j  , 2*k+1); // middle layer
    elem.v[13] = pointidx27(2*i+2, 2*j  , 2*k+1);
    elem.v[14] = pointidx27(2*i+2, 2*j+2, 2*k+1);
    elem.v[15] = pointidx27(2*i  , 2*j+2, 2*k+1);

    elem.v[16] = pointidx27(2*i+1, 2*j+0, 2*k+0); // mid segment top layer
    elem.v[17] = pointidx27(2*i+2, 2*j+1, 2*k+0);
    elem.v[18] = pointidx27(2*i+1, 2*j+2, 2*k+0);
    elem.v[19] = pointidx27(2*i+0, 2*j+1, 2*k+0);

    elem.v[20] = pointidx27(2*i+1, 2*j+1, 2*k+2);

    elem.v[21] = pointidx27(2*i+1, 2*j+0, 2*k+1); // mid segment middle layer
    elem.v[22] = pointidx27(2*i+2, 2*j+1, 2*k+1);
    elem.v[23] = pointidx27(2*i+1, 2*j+2, 2*k+1);
    elem.v[24] = pointidx27(2*i+0, 2*j+1, 2*k+1);

    elem.v[25] = pointidx27(2*i+1, 2*j+1, 2*k  );
    elem.v[26] = pointidx27(2*i+1, 2*j+1, 2*k+1);
}


void RectMesh::init_rectangular_mesh(Mesh3D& mesh)
{
    bool extvol,intvol;
    int pos_mrrs;
    Surface* smirror;

    assert(xmin<xmax);
    assert(ymin<ymax);

    apply_pml_borders();

    mesh.set_control_nodes(elem_shape);

    nelemx = int( (xmax-xmin)/xstep );
    nelemy = int( (ymax-ymin)/ystep );
    nelemz = 0;
    for(int k=0;k<nlayers;++k) {
        nelemz += nsteps[k];
    }

    printf("Creating grid mesh %d x %d x %d", nelemx, nelemy, nelemz );
    int nnodes = 0;
    if (elem_shape==8) {
        printf(" with linear elements\n");
        if (has_sph) {
            nnodes = create_linear_grid_nodes_sph(mesh);
        } else {
            nnodes = create_linear_grid_nodes(mesh);
        }
    } else if (elem_shape==27) {
        printf(" with quadratic elements\n");
        nnodes = create_quadratic_grid_nodes(mesh);
    }

    // Elements
    Elem elem(elem_shape);
    int k=0;
    if (has_mirrors) smirror = mesh.get_surface("mirror");

    for(int nl=0;nl<nlayers;++nl) {
        for(int kl=0;kl<nsteps[nl];++kl) {
            for(int j=0;j<nelemy;++j) {
                for(int i=0;i<nelemx;++i) {
                    if (elem_shape==8) {
                        create_linear_element(elem, i,j,k);
                    } else if (elem_shape==27) {
                        create_quadratic_element(elem, i,j,k);
                    }
                    for(int in=0;in<elem_shape;++in) {
                        assert(elem.v[in]<nnodes);
                        assert(elem.v[in]>=0);
                    }
                    bool W = (i<npml);
                    bool E = (i>(nelemx-npml-1));
                    bool S = (j<npml);
                    bool N = (j>(nelemy-npml-1));
                    bool U = (k<npml && nl==0);
                    bool D = (k>(nelemz-npml-1) && nl==(nlayers-1));
                    int mat = get_mat(mesh, nl, W,E,S,N,U,D);
                    // Mirrors
                    if (has_mirrors) {
                        extvol = (i<mrrs.nW || i>(nelemx-mrrs.nE-1) ||
                                       j<mrrs.nS || j>(nelemy-mrrs.nN-1) ||
                                       k<mrrs.nU || k>(nelemz-mrrs.nD-1));
                        intvol = (i>mrrs.nW && i<(nelemx-mrrs.nE-1) &&
                                       j>mrrs.nS && j<(nelemy-mrrs.nN-1) &&
                                       k>mrrs.nU && k<(nelemz-mrrs.nD-1));
                        pos_mrrs = 0;
                        if (extvol) pos_mrrs = 2;
                        if (intvol) pos_mrrs = 1;
                    }
                    // Emit element
                    if (has_mirrors) { //XXX
                        mesh.add_elem_mrrs(mat, pos_mrrs, elem);
                    } else {
                        mesh.add_elem(mat, elem);
                    }
                    // Check for free fluid surface and free pml surface
                    int dom = mesh.m_materials[mat].domain();
                    // Mirrors
                    if (has_mirrors) {
                        W = (!extvol && mrrs.W && i==mrrs.nW);
                        S = (!extvol && mrrs.S && j==mrrs.nS);
                        U = (!extvol && mrrs.U && k==mrrs.nU);
                        E = (!extvol && mrrs.E && i==(nelemx-mrrs.nE-1));
                        N = (!extvol && mrrs.N && j==(nelemy-mrrs.nN-1));
                        D = (!extvol && mrrs.D && k==(nelemz-mrrs.nD-1));
                        emit_mirr_face(smirror,dom,elem,W,E,S,N,U,D);
                    }
                }
            }
            k++;
        }
    }
}

void RectMesh::emit_free_face(Surface* surf, int dom, const Elem& elem,
                              bool W, bool E, bool S, bool N, bool U, bool D)
{
    if (W) emit_free_face(surf, dom, elem, 4);
    if (E) emit_free_face(surf, dom, elem, 2);
    if (S) emit_free_face(surf, dom, elem, 1);
    if (N) emit_free_face(surf, dom, elem, 3);
    if (U) emit_free_face(surf, dom, elem, 5);
    if (D) emit_free_face(surf, dom, elem, 0);
}

void RectMesh::emit_free_face(Surface* surf, int dom, const Elem& elem, int facenum)
{
    index_t n[4];
    for(int k=0;k<4;++k) {
        n[k] = elem.v[RefFace[facenum].v[k]];
    }
    PFace fc(n, dom);
    surf->add_face(fc, 0);

}

void RectMesh::emit_mirr_face(Surface* surf, int dom, const Elem& elem,
                              bool W, bool E, bool S, bool N, bool U, bool D)
{
    if (W) emit_mirr_face(surf, dom, elem, 4, true);
    if (E) emit_mirr_face(surf, dom, elem, 2, true);
    if (S) emit_mirr_face(surf, dom, elem, 1, false);
    if (N) emit_mirr_face(surf, dom, elem, 3, false);
    if (U) emit_mirr_face(surf, dom, elem, 5, false);
    if (D) emit_mirr_face(surf, dom, elem, 0, false);
}

void RectMesh::emit_mirr_face(Surface* surf, int dom, const Elem& elem, int facenum, bool flp)
{
    index_t n[4];
    for(int k=0;k<4;++k) {
        n[k] = elem.v[RefFace[facenum].v[k]];
    }
    PFace fc(n, dom);
    if (flp) {fc.orient = -fc.orient;}
    surf->add_face(fc, 0);

}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
