#if defined(SEM_VEC)
#define Density_(i,j,k,e) m_Density(i,j,k,e)
#else
#define Density_(i,j,k,e) m_Density(e,i,j,k)
#endif
