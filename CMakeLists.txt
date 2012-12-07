# Ce CMakeLists.txt est le projet principal permettant de construire sem2d et sem3d
# en standalone (hors mka3d)
PROJECT(SEM)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)
ENABLE_TESTING()


add_definitions(-DMKA3D)
if(OPT_COUPLAGE)
add_definitions(-DCOUPLAGE)
endif(OPT_COUPLAGE)

set(SEMDIR ${PROJECT_SOURCE_DIR})
set(SEMCOMMON ${SEMDIR}/COMMON)

set(SEM_SRC
  ${SEMCOMMON}/semdatafiles.F90
  ${SEMCOMMON}/constants.F90
  ${SEMCOMMON}/splib.F90
  ${SEMCOMMON}/invert2.F90
  ${SEMCOMMON}/invert_3d.F90
  ${SEMCOMMON}/pol_force.F90
  ${SEMCOMMON}/pol_lagrange.F90
  ${SEMCOMMON}/pow.F90
)

set(FTM_SOURCES ${SEMCOMMON}/libmpi_fantome.f90)

add_subdirectory(SEM2D)
add_subdirectory(SEM3D)
