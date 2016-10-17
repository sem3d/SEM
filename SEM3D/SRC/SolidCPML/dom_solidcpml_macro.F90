#define solidcpml_eps 1.e-12

! gamma_ab  defined after (12c) in Ref1
! gamma_abc defined after (12c) in Ref1
#define solidcpml_a0b_a1b_a2b \
    a0b = dom%Kappa(ee,0,i,j,k,bnum)*dom%Kappa(ee,1,i,j,k,bnum)*dom%Kappa(ee,2,i,j,k,bnum); \
    g0=dom%Beta(ee,0,i,j,k,bnum)-dom%Alpha(ee,0,i,j,k,bnum); \
    g1=dom%Beta(ee,1,i,j,k,bnum)-dom%Alpha(ee,1,i,j,k,bnum); \
    g2=dom%Beta(ee,2,i,j,k,bnum)-dom%Alpha(ee,2,i,j,k,bnum); \
    a1b = a0b*(g0+g1+g2); \
    g101=dom%Beta(ee,1,i,j,k,bnum)-dom%Alpha(ee,0,i,j,k,bnum)-dom%Alpha(ee,1,i,j,k,bnum); \
    g212=dom%Beta(ee,2,i,j,k,bnum)-dom%Alpha(ee,1,i,j,k,bnum)-dom%Alpha(ee,2,i,j,k,bnum); \
    g002=dom%Beta(ee,0,i,j,k,bnum)-dom%Alpha(ee,0,i,j,k,bnum)-dom%Alpha(ee,2,i,j,k,bnum); \
    a2b = a0b*(g0*g101+g1*g212+g2*g002);

#define solidcpml_a3b_a4b_a5b \
    a3b = 0.; \
    a4b = 0.; \
    a5b = 0.;
