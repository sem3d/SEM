#ifndef SOLID_FORCES_H
#define SOLID_FORCES_H


#if OPENACC
#define tFox dom%Fox(:,:,:,:,bnum)
#define tFoy dom%Foy(:,:,:,:,bnum)
#define tFoz dom%Foz(:,:,:,:,bnum)
#define tDepla dom%Depla(:,:,:,:,:,bnum)
#define tSigma dom%Sigma(:,:,:,:,:,bnum)
#define tFox_(e,i,j,k,bnum) dom%Fox(e,i,j,k,bnum)
#define tFoy_(e,i,j,k,bnum) dom%Foy(e,i,j,k,bnum)
#define tFoz_(e,i,j,k,bnum) dom%Foz(e,i,j,k,bnum)
#define tDepla_(e,i,j,k,c,bnum) dom%Depla(e,i,j,k,c,bnum)
#define tVeloc_(e,i,j,k,c,bnum) dom%Veloc(e,i,j,k,c,bnum)
#define tSigma_(e,i,j,k,c,bnum) dom%Sigma(e,i,j,k,c,bnum)
#else
#define tFox dom%Fox(:,:,:,:,0)
#define tFoy dom%Foy(:,:,:,:,0)
#define tFoz dom%Foz(:,:,:,:,0)
#define tDepla dom%Depla(:,:,:,:,:,0)
#define tSigma dom%Sigma(:,:,:,:,:,0)
#define tFox_(e,i,j,k,bnum)     dom%Fox(e,i,j,k,0)
#define tFoy_(e,i,j,k,bnum)     dom%Foy(e,i,j,k,0)
#define tFoz_(e,i,j,k,bnum)     dom%Foz(e,i,j,k,0)
#define tDepla_(e,i,j,k,c,bnum) dom%Depla(e,i,j,k,c,0)
#define tVeloc_(e,i,j,k,c,bnum) dom%Veloc(e,i,j,k,c,0)
#define tSigma_(e,i,j,k,c,bnum) dom%Sigma(e,i,j,k,c,0)
#endif


#endif
