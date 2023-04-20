#ifndef SOLID_FORCES_H
#define SOLID_FORCES_H


#if defined(OPENACC) || TEST_FORCE==1
#define tFox dom%Forces(:,:,:,:,0,bnum)
#define tFoy dom%Forces(:,:,:,:,1,bnum)
#define tFoz dom%Forces(:,:,:,:,2,bnum)
#define tDepla dom%Depla(:,:,:,:,:,bnum)
#define tSigma dom%Sigma(:,:,:,:,:,bnum)
#define tFox_(e,i,j,k,bnum) dom%Forces(e,i,j,k,0,bnum)
#define tFoy_(e,i,j,k,bnum) dom%Forces(e,i,j,k,1,bnum)
#define tFoz_(e,i,j,k,bnum) dom%Forces(e,i,j,k,2,bnum)
#define tDepla_(e,i,j,k,c,bnum) dom%Depla(e,i,j,k,c,bnum)
#define tVeloc_(e,i,j,k,c,bnum) dom%Veloc(e,i,j,k,c,bnum)
#define tSigma_(e,i,j,k,c,bnum) dom%Sigma(e,i,j,k,c,bnum)
#else
#define tFox dom%Forces(:,:,:,:,0,0)
#define tFoy dom%Forces(:,:,:,:,1,0)
#define tFoz dom%Forces(:,:,:,:,2,0)
#define tDepla dom%Depla(:,:,:,:,:,0)
#define tSigma dom%Sigma(:,:,:,:,:,0)
#define tFox_(e,i,j,k,bnum)     dom%Forces(e,i,j,k,0,0)
#define tFoy_(e,i,j,k,bnum)     dom%Forces(e,i,j,k,1,0)
#define tFoz_(e,i,j,k,bnum)     dom%Forces(e,i,j,k,2,0)
#define tDepla_(e,i,j,k,c,bnum) dom%Depla(e,i,j,k,c,0)
#define tVeloc_(e,i,j,k,c,bnum) dom%Veloc(e,i,j,k,c,0)
#define tSigma_(e,i,j,k,c,bnum) dom%Sigma(e,i,j,k,c,0)
#endif


#endif
