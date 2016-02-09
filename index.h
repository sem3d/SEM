#if defined(SEM_VEC)
#define mc_Density(i,j,k,e) mb_Density(i,j,k,e)
#else
#define mc_Density(i,j,k,e) mb_Density(e,i,j,k)
#endif
