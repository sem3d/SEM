#ifndef INDEX_H
#define INDEX_H

#if defined(SEM_VEC)
#define IND_IJKE(i,j,k,e)   e,i,j,k
#define IND_DIJKE(m,i,j,k,e)   e,m,i,j,k
#define IND_MNIJKE(m,n,i,j,k,e)   e,m,n,i,j,k

#define BEGIN_ELEM_SUBLOOP(e,e0,e1) DO e = e0,e1
#define END_ELEM_SUBLOOP  END DO
#else
#define IND_IJKE(i,j,k,e)   i,j,k,e
#define IND_DIJKE(m,i,j,k,e)   m,i,j,k,e
#define IND_MNIJKE(m,n,i,j,k,e)   m,n,i,j,k,e
#define BEGIN_ELEM_SUBLOOP(e,e0,e1) e=e0
#define END_ELEM_SUBLOOP
#endif

#define     Density_(i,j,k,e)     m_Density(IND_IJKE(i,j,k,e))
#define      Lambda_(i,j,k,e)      m_Lambda(IND_IJKE(i,j,k,e))
#define          Mu_(i,j,k,e)          m_Mu(IND_IJKE(i,j,k,e))
#define       Kappa_(i,j,k,e)       m_Kappa(IND_IJKE(i,j,k,e))
#define       Jacob_(i,j,k,e)       m_Jacob(IND_IJKE(i,j,k,e))
#define       Cij_(m,i,j,k,e)      m_Cij(IND_DIJKE(m,i,j,k,e))
#define InvGrad_(m,n,i,j,k,e) m_InvGrad(IND_MNIJKE(m,n,i,j,k,e))


#endif
