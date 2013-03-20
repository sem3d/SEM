PROJECT(mesher)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)
set(PRJ ${PROJECT_SOURCE_DIR})
set(ROOT ${CMAKE_SOURCE_DIR})

add_subdirectory(metis)

set(MESH_SRC
  ${ROOT}/COMMON/sem_hdf5.f90
  ${ROOT}/COMMON/pointlist.f90
  ${PRJ}/SRC/fich_cubit.f90
  ${PRJ}/SRC/fich_unv.f90
  ${PRJ}/SRC/local_mesh_properties.f90
  ${PRJ}/SRC/mesh2spec.f90
  ${PRJ}/SRC/mesher_mat.f90
  ${PRJ}/SRC/mesh_properties.f90
  ${PRJ}/SRC/neumann.f90
  ${PRJ}/SRC/partition_mesh.f90
  ${PRJ}/SRC/sets.f90
  ${PRJ}/SRC/solid_fluid.f90
)


add_executable(mesher ${MESH_SRC})
target_link_libraries(mesher metis ${HDF5_LIBRARIES})


