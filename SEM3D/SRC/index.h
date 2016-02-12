#ifndef INDEX_H
#define INDEX_H

/*
  Les macros definies ci-dessous doivent permettre de changer l'ordre
  des indices des champs selon la methode d'optimisation choisie.

  Il faut donc utiliser Density_(i,j,k,e) pour acceder au champ m_Density
  Selon l'option SEM_VEC, ce champ sera indexe par (i,j,k,e) ou (e,i,j,k)
  On pourrait aussi choisir de contracter i,j,k en une seule dimension
  et remplacer (i,j,k,e) par (i,1,1,e).
  
  Pour ce faire on introduira des macros pour ecrire DO I=... DO J=...

  Pour l'instant seule une boucle par element est prevue ainsi.
  
 */


#if defined(SEM_VEC)
#define IND_IJKE(i,j,k,e)   e,i,j,k
#define IND_DIJKE(m,i,j,k,e)   e,m,i,j,k
#define IND_MNIJKE(m,n,i,j,k,e)   e,m,n,i,j,k

#define BEGIN_ELEM_SUBLOOP(e,k,e0,e1) DO e = e0,e1;k=e-e0
#define END_ELEM_SUBLOOP  END DO
#else
#define IND_IJKE(i,j,k,e)   i,j,k,e
#define IND_DIJKE(m,i,j,k,e)   m,i,j,k,e
#define IND_MNIJKE(m,n,i,j,k,e)   m,n,i,j,k,e
#define BEGIN_ELEM_SUBLOOP(e,k,e0,e1) e=e0;k=0
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
