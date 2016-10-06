#include "index.h"

#define solidcpml_eps 1.e-12

#define solidcpml_x(i,xi,a,b,c,bnum,ee,mi) \
        xi = abs(dom%GlobCoord(i,dom%Idom_(a,b,c,bnum,ee)) - dom%sSubDomain(mi)%pml_pos(i));

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
        solidcpml_x    (i,xi,a,b,c,bnum,ee,mi); \
        solidcpml_kappa(i,xi,mi); \
        Li = kappa(i)
! TODO : add convolution term to Li

! (A.20) from Ref1
#define solidcpml_Lijk(i,j,k,a,b,c,bnum,ee,mi) \
        solidcpml_x    (i,xi,a,b,c,bnum,ee,mi); \
        solidcpml_kappa(i,xi,mi); \
        solidcpml_x    (j,xj,a,b,c,bnum,ee,mi); \
        solidcpml_kappa(j,xj,mi); \
        solidcpml_x    (k,xk,a,b,c,bnum,ee,mi); \
        solidcpml_kappa(k,xk,mi); \
        if (abs(kappa(k)) < solidcpml_eps) stop "ERROR: solidcpml_Lijk"; \
        Lijk = kappa(i)*kappa(j)/kappa(k)
! TODO : add convolution term to Lijk

