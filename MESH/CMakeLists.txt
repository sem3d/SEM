PROJECT(mesher)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)
set(PRJ ${PROJECT_SOURCE_DIR})
set(BASE ${CMAKE_SOURCE_DIR})
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

add_definitions(-DUSE_HDF5)
include_directories(${HDF5_INCLUDE_DIR_FORTRAN})

add_executable(mesher ${MESH_SRC})
target_link_libraries(mesher metis ${SEM_EXTRA_LIBS} )


install(TARGETS mesher DESTINATION bin)

add_executable(sem_part_grid
  ${SEMCOMMON}/h5helper.cpp
  ${PRJ}/CSRC/main_grid.cpp
  ${PRJ}/CSRC/mesh_h5_output.cpp
  ${PRJ}/CSRC/mesh.cpp
  )

add_executable(sem_part_h5
  ${SEMCOMMON}/h5helper.cpp
  ${PRJ}/CSRC/main_h5.cpp
  ${PRJ}/CSRC/mesh_h5_output.cpp
  ${PRJ}/CSRC/mesh.cpp
  ${PRJ}/CSRC/meshpart.cpp
  )



include_directories(BEFORE ${PRJ}/CSRC ${PRJ}/metis/include)
include_directories(${BASE}/COMMON)
target_link_libraries(sem_part_grid metis ${SEM_EXTRA_LIBS})
target_link_libraries(sem_part_h5 metis ${SEM_EXTRA_LIBS})

