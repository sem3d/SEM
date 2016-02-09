#if defined(SEM_VEC)
#define Density_(i,j,k,e) m_Density(i,j,k,e)
#define  Lambda_(i,j,k,e)  m_Lambda(i,j,k,e)
#define      Mu_(i,j,k,e)      m_Mu(i,j,k,e)
#define   Kappa_(i,j,k,e)   m_Kappa(i,j,k,e)
#define   Cij_(m,i,j,k,e)   m_Cij(m,i,j,k,e)
#else
#define Density_(i,j,k,e) m_Density(e,i,j,k)
#define  Lambda_(i,j,k,e)  m_Lambda(e,i,j,k)
#define      Mu_(i,j,k,e)      m_Mu(e,i,j,k)
#define   Kappa_(i,j,k,e)   m_Kappa(e,i,j,k)
#define   Cij_(m,i,j,k,e)     m_Cij(e,i,j,k,m)
#endif
