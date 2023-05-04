#ifndef LOOPS_H
#define LOOPS_H

#if VCHUNK>1
/* ATTENTION, ne pas utiliser ces macros pour l'instant, un bug du compilateur
   intel 15 fait crasher le preprocesseur... */
#define ACC_DATA(xx) !$acc data xx
#define ACC_DIR(xx)  !$acc xx
#define OMP_SIMD(xx) !$omp simd xx
#if defined(__INTEL_COMPILER)  /* IFORT */
#define OMP_DECLARE_SIMD(name,args) !$omp declare simd (name) args
#define ALIGNED(var,val) !dir$ ASSUME_ALIGNED var: val
#else /* GFORTRAN */
#define OMP_DECLARE_SIMD(name,args) !!!
#define ALIGNED(var,val) !!!
#endif
#else /* VCHUNK>1 */
#endif

#if VCHUNK>1
#define LOOP_VECTORIZE_ACC ACC_DIR(loop vector)
#define LOOP_VECTORIZE_SIMD OMP_SIMD(linear(ee))
#ifdef OPENACC
#define LOOP_VECTORIZE LOOP_VECTORIZE_ACC
#else
#define LOOP_VECTORIZE LOOP_VECTORIZE_SIMD
#endif
#define BEGIN_SUBELEM_LOOP(e,ee,bl)  do ee=0,VCHUNK-1; e = bl*VCHUNK+ee
#define BEGIN_SUBELEM_LOOP1(ee)  do ee=0,VCHUNK-1;
#define END_SUBELEM_LOOP()  enddo
#else
#define LOOP_VECTORIZE
#define BEGIN_SUBELEM_LOOP(e,ee,e0)  ee=0;e=e0;
#define BEGIN_SUBELEM_LOOP1(ee)  ee=0;
#define END_SUBELEM_LOOP()  ;
#endif


#endif LOOPS_H
