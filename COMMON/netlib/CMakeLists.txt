PROJECT(netlib)
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
ENABLE_LANGUAGE(Fortran)

set(PRJ ${PROJECT_SOURCE_DIR})

set(NETLIB_BLAS_SOURCES
${PRJ}/blas/dcopy.f
${PRJ}/blas/dgemm.f
${PRJ}/blas/dgemv.f
${PRJ}/blas/dger.f
${PRJ}/blas/dnrm2.f90
${PRJ}/blas/drot.f
${PRJ}/blas/dscal.f
${PRJ}/blas/dswap.f
${PRJ}/blas/dtrmm.f
${PRJ}/blas/dtrmv.f
${PRJ}/blas/dtrsm.f
${PRJ}/blas/lsame.f
${PRJ}/blas/xerbla.f
)
set(NETLIB_LAPACK_SOURCES
${PRJ}/lapack/la_constants.f90
${PRJ}/lapack/la_xisnan.F90
${PRJ}/lapack/dbdsqr.f
${PRJ}/lapack/dcombssq.f
${PRJ}/lapack/dgebd2.f
${PRJ}/lapack/dgebrd.f
${PRJ}/lapack/dgelq2.f
${PRJ}/lapack/dgelqf.f
${PRJ}/lapack/dgels.f
${PRJ}/lapack/dgeqr2.f
${PRJ}/lapack/dgeqrf.f
${PRJ}/lapack/dgesvd.f
${PRJ}/lapack/disnan.f
${PRJ}/lapack/dlabad.f
${PRJ}/lapack/dlabrd.f
${PRJ}/lapack/dlacpy.f
${PRJ}/lapack/dlaisnan.f
${PRJ}/lapack/dlamch.f
${PRJ}/lapack/dlange.f
${PRJ}/lapack/dlapy2.f
${PRJ}/lapack/dlarf.f
${PRJ}/lapack/dlarfb.f
${PRJ}/lapack/dlarfg.f
${PRJ}/lapack/dlarft.f
${PRJ}/lapack/dlartg.f90
${PRJ}/lapack/dlas2.f
${PRJ}/lapack/dlascl.f
${PRJ}/lapack/dlaset.f
${PRJ}/lapack/dlasq1.f
${PRJ}/lapack/dlasq2.f
${PRJ}/lapack/dlasq3.f
${PRJ}/lapack/dlasq4.f
${PRJ}/lapack/dlasq5.f
${PRJ}/lapack/dlasq6.f
${PRJ}/lapack/dlasr.f
${PRJ}/lapack/dlasrt.f
${PRJ}/lapack/dlassq.f90
${PRJ}/lapack/dlasv2.f
${PRJ}/lapack/dorg2r.f
${PRJ}/lapack/dorgbr.f
${PRJ}/lapack/dorgl2.f
${PRJ}/lapack/dorglq.f
${PRJ}/lapack/dorgqr.f
${PRJ}/lapack/dorm2r.f
${PRJ}/lapack/dormbr.f
${PRJ}/lapack/dorml2.f
${PRJ}/lapack/dormlq.f
${PRJ}/lapack/dormqr.f
${PRJ}/lapack/dtrtrs.f
${PRJ}/lapack/ieeeck.f
${PRJ}/lapack/iladlc.f
${PRJ}/lapack/iladlr.f
${PRJ}/lapack/ilaenv.f
${PRJ}/lapack/iparmq.f
${PRJ}/lapack/xerbla.f
)


add_library(netlib ${NETLIB_LAPACK_SOURCES} ${NETLIB_BLAS_SOURCES})
set_target_properties(netlib PROPERTIES
                      COMPILE_FLAGS "-I${SEMLIBDIR}/common"
                      Fortran_MODULE_DIRECTORY ${SEMLIBDIR}/common)
