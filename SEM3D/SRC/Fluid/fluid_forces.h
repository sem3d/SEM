#ifndef FLUID_FORCES_H
#define FLUID_FORCES_H


#if defined(OPENACC) || TEST_FLUID_ACC==1
#define tForcesFl  dom%ForcesFl(:,:,:,:,0,bnum)
#define tForcesFl_(e,i,j,k,bnum) dom%ForcesFl(e,i,j,k,bnum)
#define tPhi_(e,i,j,k,bnum) dom%Phi(e,i,j,k,bnum)
#else
#define tForcesFl  dom%ForcesFl(:,:,:,:,0,0)
#define tForcesFl_(e,i,j,k,bnum) dom%ForcesFl(e,i,j,k,0)
#define tPhi_(e,i,j,k,bnum) dom%Phi(e,i,j,k,0)
#endif


#endif
