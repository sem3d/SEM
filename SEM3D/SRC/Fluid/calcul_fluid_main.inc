!!! This file is part of SEM
!!!
!!! Copyright CEA, ECP, IPGP
!!!


!!! On nomme ici les fonctions communes (acc/cpu)
!!! Si une variante acc ou cpu necessite un nom de fonction differente
!!! Elle peut utiliser le meme mecanisme pour definir son nom.

#undef PROCNAME
#undef ngllvar


!!! Generic values of ngll
#ifndef NGLLVAL
#define ngllvar ngll
#define PROCNAME GLUE(PROCNAMEBASE,N)
#else
!!! Compilation time specified values of ngll
#define PROCNAME GLUE(PROCNAMEBASE,NGLLVAL)
#define ngllvar ngll0
#endif

#undef SFX
#if defined(OPENACC) || defined(OMPTARGET) || TEST_FLUID_ACC==1
#define SFX _acc
#else
#define SFX _cpu
#endif


#define SUBPROCNAME1  GLUE3(PROCNAME,_main,SFX)
#define SUBPROCNAME2  GLUE3(PROCNAME,_subloop,SFX)
#define PROC_MIRROR   GLUE(PROCNAME,_mirror)

#if defined(OPENACC) || defined(OMPTARGET) || TEST_FLUID_ACC==1
#include "calcul_forces_fluid_acc.inc"
#else
#include "calcul_forces_fluid.inc"
#endif


!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
