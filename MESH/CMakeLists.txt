PROJECT(mesher)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)
set(PRJ ${PROJECT_SOURCE_DIR})

add_subdirectory(metis)

set(MESH_SRC
  ${SEMCOMMON}/sem_hdf5.f90
  ${SEMCOMMON}/pointlist.f90
  ${PRJ}/SRC/mesh_earthchunk.f90
  ${PRJ}/SRC/fich_cubit.f90
  ${PRJ}/SRC/fich_unv.f90
  ${PRJ}/SRC/free_surf.f90
  ${PRJ}/SRC/local_mesh_properties.f90
  ${PRJ}/SRC/mesh2spec.f90
  ${PRJ}/SRC/mesher_mat.f90
  ${PRJ}/SRC/mesh_properties.f90
  ${PRJ}/SRC/neumann.f90
  ${PRJ}/SRC/partition_mesh.f90
  ${PRJ}/SRC/sets.f90
  ${PRJ}/SRC/solid_fluid.f90
  ${PRJ}/SRC/free_surf.f90
)

if(OPT_MPI)
#set(SEM_EXTRA_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
# TODO: avec mpich MPI_LIBRARY contient libmpichcxx, qui a besoin des librairies de support
# c++ et ne link pas semx.exe
set(SEM_EXTRA_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
else(OPT_MPI)
endif(OPT_MPI)

add_definitions(-DUSE_HDF5)
include_directories(${HDF5_INCLUDE_DIR_FORTRAN})

add_executable(mesher ${MESH_SRC})
target_link_libraries(mesher metis ${HDF5_LIBRARIES} )


install(TARGETS mesher DESTINATION bin)

add_executable(cmesh
  ${PRJ}/CSRC/main.cpp
  ${PRJ}/CSRC/mesh.cpp
)
include_directories(${PRJ}/CSRC ${PRJ}/metis/include)
target_link_libraries(cmesh metis ${HDF5_LIBRARIES})

