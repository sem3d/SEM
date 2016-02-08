#if defined(SEM_VEC)
#define Density(i,j,k,e) Density_(i,j,k,e)
#else
#define Density(i,j,k,e) Density_(e,i,j,k)
#endif
