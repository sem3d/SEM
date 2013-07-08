PROJECT(mesher)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)
set(PRJ ${PROJECT_SOURCE_DIR})

add_subdirectory(metis)

set(MESH_SRC
  ${SEMCOMMON}/sem_hdf5.f90
  ${SEMCOMMON}/pointlist.f90
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

set(MESH_EXTRA_LIBS "")

if (OPT_HDF5)
    add_definitions(-DUSE_HDF5)
    include_directories(${HDF5_INCLUDE_DIR_FORTRAN})
    set(MESH_EXTRA_LIBS ${MESH_EXTRA_LIBS} hdf5_fortran)
endif (OPT_HDF5)

add_executable(mesher ${MESH_SRC})
target_link_libraries(mesher metis ${MESH_EXTRA_LIBS})


install(TARGETS mesher DESTINATION bin)