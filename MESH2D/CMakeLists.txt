PROJECT(mesh2d)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)
set(PRJ ${PROJECT_SOURCE_DIR})
set(MESH_SRC
  ${SEMCOMMON}/sem_hdf5.f90
  ${PRJ}/SRC/create_2D.f90
  ${PRJ}/SRC/fich_cubit.f90 
  ${PRJ}/SRC/gid2spec_seq.f90
)

set(MESH_EXTRA_LIBS "")

if(OPT_MPI)
set(SEM_EXTRA_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
else(OPT_MPI)
endif(OPT_MPI)

add_definitions(-DUSE_HDF5)
include_directories(${HDF5_INCLUDE_DIR_FORTRAN})
set(MESH_EXTRA_LIBS ${MESH_EXTRA_LIBS} hdf5_fortran)

add_executable(mesh2d ${MESH_SRC})
target_link_libraries(mesh2d ${MESH_EXTRA_LIBS})

install(TARGETS mesh2d DESTINATION bin)
