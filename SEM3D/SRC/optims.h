#ifndef OPTIMS_H
#define OPTIMS_H

/* COMMON DECLARATIONS OPENACC OPENMP/TARGET */

#if defined(OPENMP)

#define OMPACC_ATOMIC_OP(value)  !$omp atomic value
#define OMPACC_DECL_ROUTINE_WORKER  !$omp declare target
#define OMPACC_ENTER_DATA !$omp target enter data
#define OMPACC_ENTER_DATA !$omp target exit data
#define OMPACC_COPYIN(...) !$omp&  map(to:__VA_ARGS__)
#define OMPACC_PRESENT(...) !$omp&  map(present,alloc:__VA_ARGS__)
#define OMPACC_DELETE(...) !$omp&  map(delete:__VA_ARGS__)
#define OMPACC_CREATE(...) !$omp&  map(delete:__VA_ARGS__)

#elif defined(OPENACC)

#define OMPACC_ATOMIC_OP(value)  !$acc atomic value
#define OMPACC_DECL_ROUTINE_WORKER  !$acc routine worker
#define OMPACC_ENTER_DATA !$acc enter data
#define OMPACC_EXIT_DATA !$acc enter data
#define OMPACC_COPYIN(...) !$acc&  copyin(__VA_ARGS__)
#define OMPACC_PRESENT(...) !$acc&  present(__VA_ARGS__)
#define OMPACC_DELETE(...) !$acc&  delete(__VA_ARGS__)
#define OMPACC_CREATE(...) !$acc&  create(__VA_ARGS__)

#else

#define OMPACC_ATOMIC_OP(value)
#define OMPACC_DECL_ROUTINE_WORKER

#endif



#endif // OPTIMS_H
