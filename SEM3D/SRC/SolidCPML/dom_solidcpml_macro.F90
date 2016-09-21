#include "index.h"

#define solidcpml_eps 1.e-12

#define solidcpml_x(i,xi,a,b,c,bnum,ee,mi) \
        xi = abs(dom%GlobCoord(i,dom%Idom_(a,b,c,bnum,ee)) - dom%sSubDomain(mi)%pml_pos(i));

#define solidcpml_xoverl(xyz,xi,mi) \
        xoverl = 0.; \
        if (abs(dom%sSubDomain(mi)%pml_width(xyz)) > solidcpml_eps) then; \
            xoverl = xi/abs(dom%sSubDomain(mi)%pml_width(xyz)); \
        endif;

! kappa*: (77) from Ref1
#define solidcpml_kappa(xyz,xi,mi) \
        solidcpml_xoverl(xyz,xi,mi); \
        kappa(xyz) = dom%kappa_0 + dom%kappa_1 * xoverl;

! alpha*: (76) from Ref1, beta*: (11) from Ref1, dxi: (74) from Ref1, d0: (75) from Ref1
#define solidcpml_abk(xyz,i,j,k,bnum,ee,mi) \
        solidcpml_x     (xyz,xi,i,j,k,bnum,ee,mi); \
        solidcpml_xoverl(xyz,xi,mi); \
        alpha(xyz) = dom%alphamax*(1. - xoverl); \
        solidcpml_kappa(xyz,xi,mi); \
        d0 = 0.; \
        if (abs(dom%sSubDomain(mi)%pml_width(xyz)) > solidcpml_eps) then; \
            d0 = -1.*(dom%n(xyz)+1)*dom%sSubDomain(mi)%Pspeed*log(dom%r_c); \
            d0 = d0/(2*dom%sSubDomain(mi)%pml_width(xyz)); \
        end if; \
        dxi = 0.; \
        if (abs(dom%sSubDomain(mi)%pml_width(xyz)) > solidcpml_eps) then; \
            dxi = dom%c(xyz)*d0*(xi/dom%sSubDomain(mi)%pml_width(xyz))**dom%n(xyz); \
        end if; \
        beta(xyz) = alpha(xyz) + dxi / kappa(xyz);

! gamma_abc defined after (12c) in Ref1
#define solidcpml_gamma_abc(g,a_name,a_index,b_name,b_index,c_name,c_index) \
        g = a_name(a_index) - b_name(b_index) - c_name(c_index);

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
        Lijk = kappa(i)*kappa(j)/kappa(k)
! TODO : add convolution term to Lijk

