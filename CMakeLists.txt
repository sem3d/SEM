# Ce CMakeLists.txt est le projet principal permettant de construire sem2d et sem3d
# en standalone (hors mka3d)
PROJECT(SEM)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
ENABLE_LANGUAGE(C CXX Fortran)
ENABLE_TESTING()

set(default_build_type "RelWithDebInfo")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(EXTRA_IFORT_FLAGS "")
set(EXTRA_GFORT_FLAGS "")
set(EXTRA_NVFORT_FLAGS "")

option(OPT_SINGLE "Single precision" OFF)
if (OPT_SINGLE)
  add_definitions(-DSINGLEPRECISION)
endif(OPT_SINGLE)

option(OPT_CPML "CPML" OFF)
if (OPT_CPML)
  add_definitions(-DCPML)
endif(OPT_CPML)


set(OPT_VEC "1" CACHE STRING "Vectorization chunks")

add_definitions(-DVCHUNK=${OPT_VEC})
if (${OPT_VEC}>1)
  set(EXTRA_IFORT_FLAGS  "${EXTRA_IFORT_FLAGS}  -qopenmp-simd")
  set(EXTRA_GFORT_FLAGS  "${EXTRA_GFORT_FLAGS}  -fopenmp")
  set(EXTRA_NVFORT_FLAGS "${EXTRA_NVFORT_FLAGS} -fopenmp")
endif(${OPT_VEC}>1)

option(OPT_OPENACC "OPENACC" OFF)
if (OPT_OPENACC)
  add_definitions(-DOPENACC)
  set(EXTRA_IFORT_FLAGS  "${EXTRA_IFORT_FLAGS}  ")
  set(EXTRA_GFORT_FLAGS  "${EXTRA_GFORT_FLAGS}  -fopenacc")
  set(EXTRA_NVFORT_FLAGS "${EXTRA_NVFORT_FLAGS} -acc=gpu ") # -acc=gpu -acc=multicore,... -cuda ... ?
endif(OPT_OPENACC)

option(OPT_OPENMP "OPENMP" OFF)
if (OPT_OPENMP)
  add_definitions(-DOPENMP)
  set(EXTRA_IFORT_FLAGS  "${EXTRA_IFORT_FLAGS}  -qopenmp")
  set(EXTRA_GFORT_FLAGS  "${EXTRA_GFORT_FLAGS}  -fopenmp")
  set(EXTRA_NVFORT_FLAGS "${EXTRA_NVFORT_FLAGS} ") # -acc=gpu -acc=multicore,...
endif(OPT_OPENMP)

set(OPT_MIN_GLL_OPT "4" CACHE STRING "Minimum hardcoded gll order")
set(OPT_MAX_GLL_OPT "7" CACHE STRING "Maximum hardcoded gll order")
set(OPT_MAX_GLL     "8" CACHE STRING  "Maximum gll order (any if > MAX_GLL_OPT)")

option(OPT_MPI "Compilation en mode MPI" ON)

# Obligatoire mais le flag n'est pas supprime partout
if(HDF5_DIR)
set(CMAKE_MODULE_PATH ${HDF5_DIR} ${CMAKE_MODULE_PATH})
endif(HDF5_DIR)

set(SEMDIR ${PROJECT_SOURCE_DIR})

set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 REQUIRED C Fortran)
add_definitions(-DUSE_HDF5)
include_directories(${HDF5_FORTRAN_INCLUDE_DIR} ${HDF5_INCLUDE_DIR})
##include_directories(${HDF5_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR_FORTRAN} ${HDF5_INCLUDE_DIR_CPP})
include_directories(${HDF5_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR_FORTRAN} ${HDF5_INCLUDE_DIR_CPP})
if (${HDF5_ROOT})
  include_directories(${HDF5_ROOT}/include/shared)
endif ()
include_directories(${PROJECT_BINARY_DIR}/libs/common)

if (OPT_MPI OR HDF5_IS_PARALLEL)
    set(MPI_COMPILER mpif90)
    find_package(MPI)
    include_directories(${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
    set(GLOBAL_LINK_FLAGS ${MPI_C_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS} ${MPI_Fortran_FLAGS} ${GLOBAL_LINK_FLAGS})
endif (OPT_MPI OR HDF5_IS_PARALLEL)


if (OPT_MPI)
  # Separate in case one wants to build non-mpi version of SEM against MPI enabled HDF5 library
  add_definitions(-D__MPI)
endif (OPT_MPI)

# We use c99 and c++11
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

message("Fortran COMPILER is ${CMAKE_Fortran_COMPILER_ID}")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # FLAGS ifort :
  # -fpp : preprocesseur
  # -allow nofpp_comments : indique au preprocesseur de ne pas traiter specialement les commentaires fortran
  #                         ceci permet de caser des directives openmp dans des macros
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -mkl=sequential -allow nofpp_comments ${EXTRA_IFORT_FLAGS}")
  set(USE_MKL TRUE)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -qmkl=sequential -allow nofpp_comments ${EXTRA_IFORT_FLAGS}")
  set(USE_MKL TRUE)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_NVFORT_FLAGS} -Mpreprocess ")
  set(USE_MKL FALSE)
  set(CMAKE_Fortran_DEPENDS_USE_COMPILER FALSE)
else ()
  # Assume gfortran
  # -fdefault-real-8 -fdefault-double-8
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_GFORT_FLAGS} -cpp -Wall -Wno-unused-dummy-argument -ffree-line-length-0")
  set(USE_MKL FALSE)
endif ()

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
  ${SEMCOMMON}/mat_symposdef_solver.F90 
  ${SEMCOMMON}/pow.F90
  ${SEMCOMMON}/pointlist.f90
  ${SEMCOMMON}/sem_hdf5.f90
  ${SEMCOMMON}/blas.F90
  ${SEMCOMMON}/leastsq.f90
  ${SEMCOMMON}/intset.F90
  ${SEMCOMMON}/invariants.f90
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
find_library(CEASHER cea_sherpa PATHS /usr/local/sr/lib)

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
set(SEM_CPP_LIBS ${HDF5_C_LIBRARIES})

# NETLIB/MKL
set(SEM_EXTRA_LIBS ${SEM_EXTRA_LIBS} netlib)

# HDF5
set(SEM_EXTRA_LIBS ${SEM_EXTRA_LIBS} ${HDF5_Fortran_LIBRARIES})

# MPI
if(OPT_MPI OR HDF5_IS_PARALLEL)
# TODO: avec mpich MPI_LIBRARY contient libmpichcxx, qui a besoin des librairies de support
# c++ et ne link pas semx.exe
set(SEM_EXTRA_LIBS ${SEM_EXTRA_LIBS} ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})
set(SEM_CPP_LIBS ${SEM_CPP_LIBS} ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
else(OPT_MPI OR HDF5_IS_PARALLEL)
set(SEM_SRC ${SEM_SRC} ${FTM_SOURCES})
endif(OPT_MPI OR HDF5_IS_PARALLEL)


if (CEAUSER)
add_definitions(-DHAVE_TREMAIN)
set(SEM_EXTRA_LIBS ${SEM_EXTRA_LIBS} ${CEAUSER} ${CEASHER})
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

