#include "index.h"

#define solidcpml_eps 1.e-12


! xoverl = xi/pml_width
#define solidcpml_xoverl(xyz,xi,mi) \
        xoverl = 0.; \
        if (abs(dom%sSubDomain(mi)%pml_width(xyz)) > solidcpml_eps) then; \
            xoverl = xi/abs(dom%sSubDomain(mi)%pml_width(xyz)); \
        endif; \
        if (xoverl > 1) stop "ERROR: solidcpml_xoverl";

! kappa*: (77) from Ref1
#define solidcpml_kappa(xyz,xi,mi) \
        solidcpml_xoverl(xyz,xi,mi); \
        kappa(xyz) = dom%kappa_0 + dom%kappa_1 * xoverl;

! (A.18) from Ref1
#define solidcpml_Li(i,a,b,c,bnum,ee,mi) \
        xi = abs(dom%GlobCoord(i,dom%Idom_(a,b,c,bnum,ee)) - dom%sSubDomain(mi)%pml_pos(i)); \
        solidcpml_xoverl(i,xi,mi); \
        kappa(i) = dom%kappa_0 + dom%kappa_1 * xoverl; \
        Li = kappa(i)
! TODO : add convolution term to Li

! (A.20) from Ref1
#define solidcpml_Lijk(i,j,k,a,b,c,bnum,ee,mi) \
        xi = abs(dom%GlobCoord(i,dom%Idom_(a,b,c,bnum,ee)) - dom%sSubDomain(mi)%pml_pos(i)); \
        solidcpml_xoverl(i,xi,mi); \
        kappa(i) = dom%kappa_0 + dom%kappa_1 * xoverl; \
        xj = abs(dom%GlobCoord(j,dom%Idom_(a,b,c,bnum,ee)) - dom%sSubDomain(mi)%pml_pos(j)); \
        solidcpml_xoverl(j,xj,mi); \
        kappa(j) = dom%kappa_0 + dom%kappa_1 * xoverl; \
        xk = abs(dom%GlobCoord(k,dom%Idom_(a,b,c,bnum,ee)) - dom%sSubDomain(mi)%pml_pos(k)); \
        solidcpml_xoverl(k,xk,mi); \
        kappa(k) = dom%kappa_0 + dom%kappa_1 * xoverl; \
        if (abs(kappa(k)) < solidcpml_eps) stop "ERROR: solidcpml_Lijk"; \
        Lijk = kappa(i)*kappa(j)/kappa(k)
! TODO : add convolution term to Lijk

