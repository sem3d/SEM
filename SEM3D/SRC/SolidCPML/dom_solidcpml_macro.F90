#include "index.h"

#define solidcpml_eps 1.e-12

#define solidcpml_x(i,xi,a,b,c,bnum,ee) \
        xi = abs(dom%GlobCoord(i,dom%Idom_(a,b,c,bnum,ee)) - dom%sSubDomain(bnum*VCHUNK+ee)%pml_pos(i));

#define solidcpml_xoverl(xyz,xi) \
        xoverl = 0.; \
        if (abs(dom%sSubDomain(bnum*VCHUNK+ee)%pml_width(xyz)) > solidcpml_eps) then; \
            xoverl = xi/abs(dom%sSubDomain(bnum*VCHUNK+ee)%pml_width(xyz)); \
        endif;

! kappa*: (77) from Ref1
#define solidcpml_kappa(xyz,xi) \
        solidcpml_xoverl(xyz,xi); \
        kappa(xyz) = dom%kappa_0 + dom%kappa_1 * xoverl;

! alpha*: (76) from Ref1, beta*: (11) from Ref1, dxi: (74) from Ref1, d0: (75) from Ref1
#define solidcpml_abk(xyz,i,j,k,bnum,ee) \
        solidcpml_x     (xyz,xi,i,j,k,bnum,ee); \
        solidcpml_xoverl(xyz,xi); \
        alpha(xyz) = dom%alphamax*(1. - xoverl); \
        solidcpml_kappa(xyz,xi); \
        d0 = 0.; \
        if (abs(dom%sSubDomain(bnum*VCHUNK+ee)%pml_width(xyz)) > solidcpml_eps) then; \
            d0 = -1.*(dom%n(xyz)+1)*dom%sSubDomain(bnum*VCHUNK+ee)%Pspeed*log(dom%r_c); \
            d0 = d0/(2*dom%sSubDomain(bnum*VCHUNK+ee)%pml_width(xyz)); \
        end if; \
        dxi = 0.; \
        if (abs(dom%sSubDomain(bnum*VCHUNK+ee)%pml_width(xyz)) > solidcpml_eps) then; \
            dxi = dom%c(xyz)*d0*(xi/dom%sSubDomain(bnum*VCHUNK+ee)%pml_width(xyz))**dom%n(xyz); \
        end if; \
        beta(xyz) = alpha(xyz) + dxi / kappa(xyz);

! gamma_ab defined after (12c) in Ref1
#define solidcpml_gamma_ab(g,a_name,a_index,b_name,b_index) \
        g = a_name(a_index) - b_name(b_index);

! gamma_abc defined after (12c) in Ref1
#define solidcpml_gamma_abc(g,a_name,a_index,b_name,b_index,c_name,c_index) \
        g = a_name(a_index) - b_name(b_index) - c_name(c_index);

! (A.18) from Ref1
#define solidcpml_Li(i,a,b,c,bnum,ee) \
        solidcpml_x    (i,xi,a,b,c,bnum,ee); \
        solidcpml_kappa(i,xi); \
        Li = kappa(i)
! TODO : add convolution term to Li

! (A.20) from Ref1
#define solidcpml_Lijk(i,j,k,a,b,c,bnum,ee) \
        solidcpml_x    (i,xi,a,b,c,bnum,ee); \
        solidcpml_kappa(i,xi); \
        solidcpml_x    (j,xj,a,b,c,bnum,ee); \
        solidcpml_kappa(j,xj); \
        solidcpml_x    (k,xk,a,b,c,bnum,ee); \
        solidcpml_kappa(k,xk); \
        Lijk = kappa(i)*kappa(j)/kappa(k)
! TODO : add convolution term to Lijk

