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
/* ATTENTION, ne pas utiliser ces macros pour l'instant, un bug du compilateur
   intel 15 fait crasher le preprocesseur... */
#define OMP_SIMD(xx) !$omp simd xx
#if defined(__INTEL_COMPILER)  /* IFORT */
#define OMP_DECLARE_SIMD(name,args) !$omp declare simd (name) args
#define ALIGNED(var,val) !dir$ ASSUME_ALIGNED var: val
#else /* GFORTRAN */
#define OMP_DECLARE_SIMD(name,args) !!!
#define ALIGNED(var,val) !!!
#endif
#else /* SEM_VEC */
#endif
#if defined(SEM_VEC)
#define CHUNK  64
#define IND_IJKE(i,j,k,e)           e,i,j,k
#define IND_MNE(m,n,e)           e,m,n
#define IND_DIJKE(m,i,j,k,e)      e,m,i,j,k
#define IND_MNIJKE(m,n,i,j,k,e) e,m,n,i,j,k
#define IND_IJKNE(i,j,k,n,e)      e,n,i,j,k
#define IND_IJKE(i,j,k,e)           e,i,j,k
#define IND_NIJKE(n,i,j,k,e)        e,i,j,k,n
#define SUBELEM_LOOP_DIR  !dir$ simd
#define BEGIN_SUBELEM_LOOP(e,ee,e0)  do ee=0,CHUNK-1; e = e0+ee

#define END_SUBELEM_LOOP()  enddo
#else
#define CHUNK 1
#define IND_IJKE(i,j,k,e)           i,j,k,e
#define IND_MNE(m,n,e)           m,n,e
#define IND_DIJKE(m,i,j,k,e)      m,i,j,k,e
#define IND_MNIJKE(m,n,i,j,k,e) m,n,i,j,k,e
#define IND_IJKNE(i,j,k,n,e)      n,i,j,k,e
#define IND_IJKE(i,j,k,e)           i,j,k,e
#define IND_NIJKE(n,i,j,k,e)        i,j,k,n,e
#define SUBELEM_LOOP_DIR
#define BEGIN_SUBELEM_LOOP(e,ee,e0)  ee=0;e=e0;
#define END_SUBELEM_LOOP()  ;
#endif

#define     Density_(i,j,k,e)     m_Density(IND_IJKE(i,j,k,e))
#define    IDensity_(i,j,k,e)    m_IDensity(IND_IJKE(i,j,k,e))
#define      Lambda_(i,j,k,e)      m_Lambda(IND_IJKE(i,j,k,e))
#define          Mu_(i,j,k,e)          m_Mu(IND_IJKE(i,j,k,e))
#define       Kappa_(i,j,k,e)       m_Kappa(IND_IJKE(i,j,k,e))
#define       Jacob_(i,j,k,e)       m_Jacob(IND_IJKE(i,j,k,e))
#define        Idom_(i,j,k,e)        m_Idom(IND_IJKE(i,j,k,e))

#define       Cij_(m,i,j,k,e)      m_Cij(IND_DIJKE(m,i,j,k,e))
#define InvGrad_(m,n,i,j,k,e) m_InvGrad(IND_MNIJKE(m,n,i,j,k,e))

#define  PMLVeloc_(i,j,k,n,e)  m_PMLVeloc(IND_IJKNE(i,j,k,n,e))
#define PMLDumpSx_(i,j,k,n,e) m_PMLDumpSx(IND_IJKNE(i,j,k,n,e))
#define PMLDumpSy_(i,j,k,n,e) m_PMLDumpSy(IND_IJKNE(i,j,k,n,e))
#define PMLDumpSz_(i,j,k,n,e) m_PMLDumpSz(IND_IJKNE(i,j,k,n,e))

#define                 Q_(i,j,k,e)                  m_Q(IND_IJKE(i,j,k,e))
#define                Qs_(i,j,k,e)                 m_Qs(IND_IJKE(i,j,k,e))
#define                Qp_(i,j,k,e)                 m_Qp(IND_IJKE(i,j,k,e))
#define        epsilonvol_(i,j,k,e)         m_epsilonvol(IND_IJKE(i,j,k,e))
#define     epsilondev_xx_(i,j,k,e)      m_epsilondev_xx(IND_IJKE(i,j,k,e))
#define     epsilondev_yy_(i,j,k,e)      m_epsilondev_yy(IND_IJKE(i,j,k,e))
#define     epsilondev_xy_(i,j,k,e)      m_epsilondev_xy(IND_IJKE(i,j,k,e))
#define     epsilondev_xz_(i,j,k,e)      m_epsilondev_xz(IND_IJKE(i,j,k,e))
#define     epsilondev_yz_(i,j,k,e)      m_epsilondev_yz(IND_IJKE(i,j,k,e))
#define            R_xx_(n,i,j,k,e)            m_R_xx(IND_NIJKE(n,i,j,k,e))
#define            R_yy_(n,i,j,k,e)            m_R_yy(IND_NIJKE(n,i,j,k,e))
#define            R_xy_(n,i,j,k,e)            m_R_xy(IND_NIJKE(n,i,j,k,e))
#define            R_xz_(n,i,j,k,e)            m_R_xz(IND_NIJKE(n,i,j,k,e))
#define            R_yz_(n,i,j,k,e)            m_R_yz(IND_NIJKE(n,i,j,k,e))
#define           R_vol_(n,i,j,k,e)           m_R_vol(IND_NIJKE(n,i,j,k,e))

#define     omega_tau_s_(n,i,j,k,e)     m_omega_tau_s(IND_NIJKE(n,i,j,k,e))
#define       agamma_mu_(n,i,j,k,e)       m_agamma_mu(IND_NIJKE(n,i,j,k,e))
#define    agamma_kappa_(n,i,j,k,e)    m_agamma_kappa(IND_NIJKE(n,i,j,k,e))

#define         onemPbeta_(i,j,k,e)          m_onemPbeta(IND_IJKE(i,j,k,e))
#define         onemSbeta_(i,j,k,e)          m_onemSbeta(IND_IJKE(i,j,k,e))
#define factor_common_3_(n,i,j,k,e) m_factor_common_3(IND_NIJKE(n,i,j,k,e))
#define      alphaval_3_(n,i,j,k,e)      m_alphaval_3(IND_NIJKE(n,i,j,k,e))
#define       betaval_3_(n,i,j,k,e)       m_betaval_3(IND_NIJKE(n,i,j,k,e))
#define      gammaval_3_(n,i,j,k,e)      m_gammaval_3(IND_NIJKE(n,i,j,k,e))
#define factor_common_P_(n,i,j,k,e) m_factor_common_P(IND_NIJKE(n,i,j,k,e))
#define      alphaval_P_(n,i,j,k,e)      m_alphaval_P(IND_NIJKE(n,i,j,k,e))
#define       betaval_P_(n,i,j,k,e)       m_betaval_P(IND_NIJKE(n,i,j,k,e))
#define      gammaval_P_(n,i,j,k,e)      m_gammaval_P(IND_NIJKE(n,i,j,k,e))

#define Diagonal_Stress1_(i,j,k,n,e) m_Diagonal_Stress1(IND_IJKNE(i,j,k,n,e))
#define Diagonal_Stress2_(i,j,k,n,e) m_Diagonal_Stress2(IND_IJKNE(i,j,k,n,e))
#define Diagonal_Stress3_(i,j,k,n,e) m_Diagonal_Stress3(IND_IJKNE(i,j,k,n,e))
#define  Diagonal_Stress_(i,j,k,n,e)  m_Diagonal_Stress(IND_IJKNE(i,j,k,n,e))
#define Residual_Stress1_(i,j,k,n,e) m_Residual_Stress1(IND_IJKNE(i,j,k,n,e))
#define Residual_Stress2_(i,j,k,n,e) m_Residual_Stress2(IND_IJKNE(i,j,k,n,e))
#define Residual_Stress3_(i,j,k,n,e) m_Residual_Stress3(IND_IJKNE(i,j,k,n,e))
#define  Residual_Stress_(i,j,k,n,e)  m_Residual_Stress(IND_IJKNE(i,j,k,n,e))
#define        PMLDumpSx_(i,j,k,n,e)        m_PMLDumpSx(IND_IJKNE(i,j,k,n,e))
#define        PMLDumpSy_(i,j,k,n,e)        m_PMLDumpSy(IND_IJKNE(i,j,k,n,e))
#define        PMLDumpSz_(i,j,k,n,e)        m_PMLDumpSz(IND_IJKNE(i,j,k,n,e))

#endif

#define part_deriv_ijke(Var,d,dS_dxi,dS_deta,dS_dzeta,dxx,dxy,dxz) \
        dS_dxi   = 0.0D+0; \
        dS_deta  = 0.0D+0; \
        dS_dzeta = 0.0D+0; \
        DO L = 0, ngll-1;  \
            dS_dxi   = dS_dxi  +Var(ee,L,J,K,d)*dom%hprime(L,I); \
            dS_deta  = dS_deta +Var(ee,I,L,K,d)*dom%hprime(L,J); \
            dS_dzeta = dS_dzeta+Var(ee,I,J,L,d)*dom%hprime(L,K); \
        END DO; \
        dxx = dS_dxi*InvGrad(IND_MNE(0,0,ee))+dS_deta*InvGrad(IND_MNE(0,1,ee))+dS_dzeta*InvGrad(IND_MNE(0,2,ee)); \
        dxy = dS_dxi*InvGrad(IND_MNE(1,0,ee))+dS_deta*InvGrad(IND_MNE(1,1,ee))+dS_dzeta*InvGrad(IND_MNE(1,2,ee)); \
        dxz = dS_dxi*InvGrad(IND_MNE(2,0,ee))+dS_deta*InvGrad(IND_MNE(2,1,ee))+dS_dzeta*InvGrad(IND_MNE(2,2,ee));


#define part2_deriv_ijke(Var,d,dS_dxi,dS_deta,dS_dzeta,dxx,dxy,dxz) \
        dS_dxi   = 0.0D+0; \
        dS_deta  = 0.0D+0; \
        dS_dzeta = 0.0D+0; \
        DO L = 0, ngll-1;  \
            dS_dxi   = dS_dxi  +Var(ee,L,J,K,d)*dom%hprime(L,I); \
            dS_deta  = dS_deta +Var(ee,I,L,K,d)*dom%hprime(L,J); \
            dS_dzeta = dS_dzeta+Var(ee,I,J,L,d)*dom%hprime(L,K); \
        END DO; \
        dxx = dS_dxi*dom%InvGrad_(0,0,i,j,k,e))+dS_deta*dom%InvGrad_(0,1,i,j,k,e))+dS_dzeta*dom%InvGrad_(0,2,i,j,k,e)); \
        dxy = dS_dxi*dom%InvGrad_(1,0,i,j,k,e))+dS_deta*dom%InvGrad_(1,1,i,j,k,e))+dS_dzeta*dom%InvGrad_(1,2,i,j,k,e)); \
        dxz = dS_dxi*dom%InvGrad_(2,0,i,j,k,e))+dS_deta*dom%InvGrad_(2,1,i,j,k,e))+dS_dzeta*dom%InvGrad_(2,2,i,j,k,e));


#define local_deriv_ijke(Var,d,dS_dxi,dS_deta,dS_dzeta) \
        dS_dxi   = 0.0D+0; \
        dS_deta  = 0.0D+0; \
        dS_dzeta = 0.0D+0; \
        DO L = 0, ngll-1;  \
            dS_dxi   = dS_dxi  +Var(ee,L,J,K,d)*dom%hprime(L,I); \
            dS_deta  = dS_deta +Var(ee,I,L,K,d)*dom%hprime(L,J); \
            dS_dzeta = dS_dzeta+Var(ee,I,J,L,d)*dom%hprime(L,K); \
        END DO;
