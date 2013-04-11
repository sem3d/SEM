# Ce CMakeLists.txt est le projet principal permettant de construire sem2d et sem3d
# en standalone (hors mka3d)
PROJECT(SEM)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(C Fortran)
ENABLE_TESTING()


if (NOT BUILD_WITH_MKA)
# Hors mka il faut chercher MPI et HDF5 ici
option(OPT_MPI "Compilation en mode MPI")

# Obligatoire mais le flag n'est pas supprime partout
find_package(HDF5 REQUIRED C Fortran)
add_definitions(-DUSE_HDF5)
include_directories(${HDF5_Fortran_INCLUDE_PATH} ${HDF5_C_INCLUDE_DIR})

if (OPT_MPI)
    set(MPI_COMPILER mpif90)
    find_package(MPI)
    add_definitions(-D__MPI)
    include_directories(${MPI_INCLUDE_PATH})
    set(GLOBAL_LINK_FLAGS ${MPI_LINK_FLAGS} ${GLOBAL_LINK_FLAGS})
endif (OPT_MPI)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -fpp -mkl=sequential")
  set(USE_MKL TRUE)
else (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Assume gfortran
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0")
  set(USE_MKL FALSE)
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
endif (NOT BUILD_WITH_MKA)


add_definitions(-DMKA3D)
if(OPT_COUPLAGE)
add_definitions(-DCOUPLAGE)
endif(OPT_COUPLAGE)

set(SEMDIR ${PROJECT_SOURCE_DIR})
set(SEMCOMMON ${SEMDIR}/COMMON)

set(SEM_SRC
  ${SEMCOMMON}/semdatafiles.F90
  ${SEMCOMMON}/protrep.f90
  ${SEMCOMMON}/sem_c_config.F90
  ${SEMCOMMON}/constants.F90
  ${SEMCOMMON}/splib.F90
  ${SEMCOMMON}/invert2.F90
  ${SEMCOMMON}/invert_3d.F90
  ${SEMCOMMON}/pol_force.F90
  ${SEMCOMMON}/lagrange_prop.f90
  ${SEMCOMMON}/pow.F90
  ${SEMCOMMON}/pointlist.f90
)

set(SEM_SRC_C
  ${SEMCOMMON}/read_input.c
  ${SEMCOMMON}/file_scan.c
  ${SEMCOMMON}/sem_input.c
  ${SEMCOMMON}/tremain_c.c
)

if(NOT BUILD_WITH_MKA)
set(SEM_SRC ${SEM_SRC}
  ${SEMCOMMON}/stub_tremain.F90
)
endif(NOT BUILD_WITH_MKA)


set(FTM_SOURCES ${SEMCOMMON}/libmpi_fantome.f90)

add_subdirectory(SEM2D)
add_subdirectory(SEM3D)
add_subdirectory(MESH)


add_executable(test_sem_path ${SEMCOMMON}/test_sem_path.f90 ${SEMCOMMON}/semdatafiles.F90)
