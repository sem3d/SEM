# Ce CMakeLists.txt est le projet principal permettant de construire sem2d et sem3d
# en standalone (hors mka3d)
PROJECT(SEM)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(C CXX Fortran)
ENABLE_TESTING()

option(OPT_CPML "CPML" OFF)
if (OPT_CPML)
  add_definitions(-DCPML)
endif(OPT_CPML)

set(OPT_VEC "1" CACHE STRING "Vectorization chunks")

add_definitions(-DVCHUNK=${OPT_VEC})
if (${OPT_VEC}>1)
  set(EXTRA_IFORT_FLAGS "-qopenmp-simd")
  set(EXTRA_GFORT_FLAGS "-fopenmp")
endif(${OPT_VEC}>1)

if (NOT BUILD_WITH_MKA)
# Hors mka il faut chercher MPI et HDF5 ici
option(OPT_MPI "Compilation en mode MPI" ON)

# Obligatoire mais le flag n'est pas supprime partout
if(HDF5_DIR)
set(CMAKE_MODULE_PATH ${HDF5_DIR} ${CMAKE_MODULE_PATH})
endif(HDF5_DIR)

set(SEMDIR ${PROJECT_SOURCE_DIR})

find_package(HDF5 REQUIRED C Fortran)
add_definitions(-DUSE_HDF5)
include_directories(${HDF5_FORTRAN_INCLUDE_DIR} ${HDF5_C_INCLUDE_DIR} ${HDF5_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR} ${HDF5_INCLUDE_DIR_FORTRAN} ${HDF5_INCLUDE_DIR_CPP})
include_directories(${PROJECT_BINARY_DIR}/libs/common)

if (OPT_MPI OR HDF5_IS_PARALLEL)
    set(MPI_COMPILER mpif90)
    find_package(MPI)
    include_directories(${MPI_INCLUDE_PATH})
    set(GLOBAL_LINK_FLAGS ${MPI_LINK_FLAGS} ${GLOBAL_LINK_FLAGS})
endif (OPT_MPI OR HDF5_IS_PARALLEL)


if (OPT_MPI)
  # Separate in case one wants to build non-mpi version of SEM against MPI enabled HDF5 library
  add_definitions(-D__MPI)
endif (OPT_MPI)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # FLAGS ifort :
  # -fpp : preprocesseur
  # -allow nofpp_comments : indique au preprocesseur de ne pas traiter specialement les commentaires fortran
  #                         ceci permet de caser des directives openmp dans des macros
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -fpp -mkl=sequential -allow nofpp_comments ${EXTRA_IFORT_FLAGS}")
  set(USE_MKL TRUE)
else (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Assume gfortran
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_GFORT_FLAGS} -cpp -Wall -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0")
  set(USE_MKL FALSE)
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
endif (NOT BUILD_WITH_MKA)

add_definitions(-DMKA3D)

message(STATUS "${CMAKE_CXX_FLAGS}")

set(SEMCOMMON ${SEMDIR}/COMMON)

set(SEM_SRC
  ${SEMCOMMON}/sem_c_bindings.F90
  ${SEMCOMMON}/semdatafiles.F90
  ${SEMCOMMON}/protrep.f90
  ${SEMCOMMON}/constants.F90
  ${SEMCOMMON}/sem_c_config.F90
  ${SEMCOMMON}/splib.F90
  ${SEMCOMMON}/invert2.F90
  ${SEMCOMMON}/invert_3d.F90
  ${SEMCOMMON}/pol_force.F90
  ${SEMCOMMON}/lagrange_prop.f90
  ${SEMCOMMON}/pow.F90
  ${SEMCOMMON}/pointlist.f90
  ${SEMCOMMON}/sem_hdf5.f90
  ${SEMCOMMON}/blas.F90
  ${SEMCOMMON}/leastsq.f90
  ${SEMCOMMON}/intset.F90
)

set(SEM_SRC_C
  ${SEMCOMMON}/read_input.c
  ${SEMCOMMON}/read_material.c
  ${SEMCOMMON}/file_scan.c
  ${SEMCOMMON}/sem_input.c
  ${SEMCOMMON}/tremain_c.c
  ${SEMCOMMON}/fileutils.c
)

find_library(CEAUSER cea_user PATHS /usr/local/sr/lib)

if(CEAUSER STREQUAL "CEAUSER-NOTFOUND")
set(SEM_SRC ${SEM_SRC}
  ${SEMCOMMON}/stub_tremain.F90
)
endif(CEAUSER STREQUAL "CEAUSER-NOTFOUND")


set(FTM_SOURCES ${SEMCOMMON}/libmpi_fantome.f90)
# This command is used to find a library. A cache entry named by <VAR> is created to store the result of this command. If the library is found the result is stored in the variable and the search will not be repeated unless the variable is cleared. If nothing is found, the re
# On veut compiler les .f90 avec le preprocesseur cpp
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set_source_files_properties(${SEM3D_SOURCES_F} ${FTM_SOURCES} PROPERTIES COMPILE_FLAGS "-fpp")
else (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set_source_files_properties(${SEM3D_SOURCES_F} ${FTM_SOURCES} PROPERTIES COMPILE_FLAGS "-x f95-cpp-input")
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")


set(SEM_EXTRA_LIBS)

# NETLIB/MKL
set(SEM_EXTRA_LIBS ${SEM_EXTRA_LIBS} netlib)

# HDF5
set(SEM_EXTRA_LIBS ${SEM_EXTRA_LIBS} ${HDF5_LIBRARIES})

# MPI
if(OPT_MPI OR HDF5_IS_PARALLEL)
# TODO: avec mpich MPI_LIBRARY contient libmpichcxx, qui a besoin des librairies de support
# c++ et ne link pas semx.exe
set(SEM_EXTRA_LIBS ${SEM_EXTRA_LIBS} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
else(OPT_MPI OR HDF5_IS_PARALLEL)
set(SEM_SRC ${SEM_SRC} ${FTM_SOURCES})
endif(OPT_MPI OR HDF5_IS_PARALLEL)


if (CEAUSER)
add_definitions(-DHAVE_TREMAIN)
set(SEM_EXTRA_LIBS ${SEM_EXTRA_LIBS} ${CEAUSER})
endif(CEAUSER)

message("SEMLIBS: '${SEM_EXTRA_LIBS}'")

set(SEMLIBDIR ${PROJECT_BINARY_DIR}/libs)
#SEM COMMON LIBRARY
add_library(semcommon ${SEM_SRC} ${SEM_SRC_C})
set_target_properties(semcommon PROPERTIES Fortran_MODULE_DIRECTORY ${SEMLIBDIR}/common)
set(SEM_EXTRA_LIBS semcommon ${SEM_EXTRA_LIBS})

#RANDOM FIELD LIBRARY
option(OPT_RF "RF" OFF)
if (OPT_RF)
  add_subdirectory(randomField)
endif (OPT_RF)

add_subdirectory(COMMON/netlib)
option(OPT_SEM2D "Compile SEM2D" OFF) # Save compilation time (in particular with intel : save half an hour !...)
if(OPT_SEM2D)
  add_subdirectory(MESH2D)
  add_subdirectory(SEM2D)
  SET_TARGET_PROPERTIES(sem2d.exe PROPERTIES LINKER_LANGUAGE Fortran)
endif(OPT_SEM2D)
add_subdirectory(SEM3D)
add_subdirectory(MESH)

add_executable(test_sem_path ${SEMCOMMON}/test_sem_path.f90   ${SEMCOMMON}/fileutils.c ${SEMCOMMON}/sem_c_bindings.F90 ${SEMCOMMON}/semdatafiles.F90)

SET_TARGET_PROPERTIES(sem3d.exe PROPERTIES LINKER_LANGUAGE Fortran)

