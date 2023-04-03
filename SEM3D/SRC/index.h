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
#if VCHUNK>1
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
#else /* VCHUNK>1 */
#endif


#define IND_IJKE(i,j,k,eb,ec)        ec,i,j,k,eb
#define IND_MNE(m,n,eb,ec)           ec,m,n,eb
#define IND_NE(n,eb,ec)              ec,n,eb
#define IND_DIJKE(m,i,j,k,eb,ec)     ec,m,i,j,k,eb
#define IND_MNIJKE(m,n,i,j,k,eb,ec)  ec,m,n,i,j,k,eb
#define IND_IJKNE(i,j,k,n,eb,ec)     ec,n,i,j,k,eb
#define IND_IJKE(i,j,k,eb,ec)        ec,i,j,k,eb
#define IND_NIJKE(n,i,j,k,eb,ec)     ec,n,i,j,k,eb

#if VCHUNK>1
#define BEGIN_SUBELEM_LOOP(e,ee,bl)  do ee=0,VCHUNK-1; e = bl*VCHUNK+ee
#define END_SUBELEM_LOOP()  enddo
#else
#define BEGIN_SUBELEM_LOOP(e,ee,e0)  ee=0;e=e0;
#define END_SUBELEM_LOOP()  ;
#endif

#define     Density_(i,j,k,eb,ec)     m_Density(IND_IJKE(i,j,k,eb,ec))
#define    IDensity_(i,j,k,eb,ec)    m_IDensity(IND_IJKE(i,j,k,eb,ec))
#define      Lambda_(i,j,k,eb,ec)      m_Lambda(IND_IJKE(i,j,k,eb,ec))
#define          Mu_(i,j,k,eb,ec)          m_Mu(IND_IJKE(i,j,k,eb,ec))
#define       Kappa_(i,j,k,eb,ec)       m_Kappa(IND_IJKE(i,j,k,eb,ec))
#define       Jacob_(i,j,k,eb,ec)       m_Jacob(IND_IJKE(i,j,k,eb,ec))
#define        Idom_(i,j,k,eb,ec)        m_Idom(IND_IJKE(i,j,k,eb,ec))
#define        Itrace_(i,j,eb,ec)        m_Itrace(IND_MNE(i,j,eb,ec))

#define        syld_(i,j,k,eb,ec)        m_syld(IND_IJKE(i,j,k,eb,ec))
#define        ckin_(i,j,k,eb,ec)        m_ckin(IND_IJKE(i,j,k,eb,ec))
#define        kkin_(i,j,k,eb,ec)        m_kkin(IND_IJKE(i,j,k,eb,ec))
#define        rinf_(i,j,k,eb,ec)        m_rinf(IND_IJKE(i,j,k,eb,ec))
#define        biso_(i,j,k,eb,ec)        m_biso(IND_IJKE(i,j,k,eb,ec))
#define      radius_(i,j,k,eb,ec)      m_radius(IND_IJKE(i,j,k,eb,ec))
#define    strain_(n,i,j,k,eb,ec)   m_strain(IND_NIJKE(n,i,j,k,eb,ec))
#define    stress_(n,i,j,k,eb,ec)   m_stress(IND_NIJKE(n,i,j,k,eb,ec))
#define    center_(n,i,j,k,eb,ec)   m_center(IND_NIJKE(n,i,j,k,eb,ec))
#define  plstrain_(n,i,j,k,eb,ec) m_plstrain(IND_NIJKE(n,i,j,k,eb,ec))

#define       Cij_(m,i,j,k,eb,ec)      m_Cij(IND_DIJKE(m,i,j,k,eb,ec))
#define InvGrad_(m,n,i,j,k,eb,ec) m_InvGrad(IND_MNIJKE(m,n,i,j,k,eb,ec))

#define  PMLVeloc_(i,j,k,n,eb,ec)  m_PMLVeloc(IND_IJKNE(i,j,k,n,eb,ec))
#define PMLDumpSx_(i,j,k,n,eb,ec) m_PMLDumpSx(IND_IJKNE(i,j,k,n,eb,ec))
#define PMLDumpSy_(i,j,k,n,eb,ec) m_PMLDumpSy(IND_IJKNE(i,j,k,n,eb,ec))
#define PMLDumpSz_(i,j,k,n,eb,ec) m_PMLDumpSz(IND_IJKNE(i,j,k,n,eb,ec))

#define                Qs_(i,j,k,eb,ec)                 m_Qs(IND_IJKE(i,j,k,eb,ec))
#define                Qp_(i,j,k,eb,ec)                 m_Qp(IND_IJKE(i,j,k,eb,ec))
#define        epsilonvol_(i,j,k,eb,ec)         m_epsilonvol(IND_IJKE(i,j,k,eb,ec))
#define     epsilondev_xx_(i,j,k,eb,ec)      m_epsilondev_xx(IND_IJKE(i,j,k,eb,ec))
#define     epsilondev_yy_(i,j,k,eb,ec)      m_epsilondev_yy(IND_IJKE(i,j,k,eb,ec))
#define     epsilondev_xy_(i,j,k,eb,ec)      m_epsilondev_xy(IND_IJKE(i,j,k,eb,ec))
#define     epsilondev_xz_(i,j,k,eb,ec)      m_epsilondev_xz(IND_IJKE(i,j,k,eb,ec))
#define     epsilondev_yz_(i,j,k,eb,ec)      m_epsilondev_yz(IND_IJKE(i,j,k,eb,ec))
#define            R_xx_(n,i,j,k,eb,ec)            m_R_xx(IND_NIJKE(n,i,j,k,eb,ec))
#define            R_yy_(n,i,j,k,eb,ec)            m_R_yy(IND_NIJKE(n,i,j,k,eb,ec))
#define            R_xy_(n,i,j,k,eb,ec)            m_R_xy(IND_NIJKE(n,i,j,k,eb,ec))
#define            R_xz_(n,i,j,k,eb,ec)            m_R_xz(IND_NIJKE(n,i,j,k,eb,ec))
#define            R_yz_(n,i,j,k,eb,ec)            m_R_yz(IND_NIJKE(n,i,j,k,eb,ec))
#define           R_vol_(n,i,j,k,eb,ec)           m_R_vol(IND_NIJKE(n,i,j,k,eb,ec))

#define     omega_tau_s_(n,i,j,k,eb,ec)     m_omega_tau_s(IND_NIJKE(n,i,j,k,eb,ec))
#define       agamma_mu_(n,i,j,k,eb,ec)       m_agamma_mu(IND_NIJKE(n,i,j,k,eb,ec))
#define    agamma_kappa_(n,i,j,k,eb,ec)    m_agamma_kappa(IND_NIJKE(n,i,j,k,eb,ec))

#define         onemPbeta_(i,j,k,eb,ec)          m_onemPbeta(IND_IJKE(i,j,k,eb,ec))
#define         onemSbeta_(i,j,k,eb,ec)          m_onemSbeta(IND_IJKE(i,j,k,eb,ec))
#define factor_common_3_(n,i,j,k,eb,ec) m_factor_common_3(IND_NIJKE(n,i,j,k,eb,ec))
#define      alphaval_3_(n,i,j,k,eb,ec)      m_alphaval_3(IND_NIJKE(n,i,j,k,eb,ec))
#define       betaval_3_(n,i,j,k,eb,ec)       m_betaval_3(IND_NIJKE(n,i,j,k,eb,ec))
#define      gammaval_3_(n,i,j,k,eb,ec)      m_gammaval_3(IND_NIJKE(n,i,j,k,eb,ec))
#define factor_common_P_(n,i,j,k,eb,ec) m_factor_common_P(IND_NIJKE(n,i,j,k,eb,ec))
#define      alphaval_P_(n,i,j,k,eb,ec)      m_alphaval_P(IND_NIJKE(n,i,j,k,eb,ec))
#define       betaval_P_(n,i,j,k,eb,ec)       m_betaval_P(IND_NIJKE(n,i,j,k,eb,ec))
#define      gammaval_P_(n,i,j,k,eb,ec)      m_gammaval_P(IND_NIJKE(n,i,j,k,eb,ec))

#define         Stress1_(i,j,k,n,eb,ec)         m_Stress1(IND_IJKNE(i,j,k,n,eb,ec))
#define         Stress2_(i,j,k,n,eb,ec)         m_Stress2(IND_IJKNE(i,j,k,n,eb,ec))
#define         Stress3_(i,j,k,n,eb,ec)         m_Stress3(IND_IJKNE(i,j,k,n,eb,ec))
#define        PMLDumpSx_(i,j,k,n,eb,ec)        m_PMLDumpSx(IND_IJKNE(i,j,k,n,eb,ec))
#define        PMLDumpSy_(i,j,k,n,eb,ec)        m_PMLDumpSy(IND_IJKNE(i,j,k,n,eb,ec))
#define        PMLDumpSz_(i,j,k,n,eb,ec)        m_PMLDumpSz(IND_IJKNE(i,j,k,n,eb,ec))

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
        dxx = dS_dxi*dom%InvGrad_(0,0,i,j,k,bnum,ee)+dS_deta*dom%InvGrad_(0,1,i,j,k,bnum,ee)+dS_dzeta*dom%InvGrad_(0,2,i,j,k,bnum,ee); \
        dxy = dS_dxi*dom%InvGrad_(1,0,i,j,k,bnum,ee)+dS_deta*dom%InvGrad_(1,1,i,j,k,bnum,ee)+dS_dzeta*dom%InvGrad_(1,2,i,j,k,bnum,ee); \
        dxz = dS_dxi*dom%InvGrad_(2,0,i,j,k,bnum,ee)+dS_deta*dom%InvGrad_(2,1,i,j,k,bnum,ee)+dS_dzeta*dom%InvGrad_(2,2,i,j,k,bnum,ee);


#define local_deriv_ijke(Var,d,dS_dxi,dS_deta,dS_dzeta) \
        dS_dxi   = 0.0D+0; \
        dS_deta  = 0.0D+0; \
        dS_dzeta = 0.0D+0; \
        DO L = 0, ngll-1;  \
            dS_dxi   = dS_dxi  +Var(ee,L,J,K,d)*dom%hprime(L,I); \
            dS_deta  = dS_deta +Var(ee,I,L,K,d)*dom%hprime(L,J); \
            dS_dzeta = dS_dzeta+Var(ee,I,J,L,d)*dom%hprime(L,K); \
        END DO;
        
