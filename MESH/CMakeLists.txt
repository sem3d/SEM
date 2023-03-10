PROJECT(mesher)
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
ENABLE_LANGUAGE(C CXX Fortran)
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

#add_executable(oldmesher ${MESH_SRC})
#target_link_libraries(oldmesher metis ${SEM_EXTRA_LIBS} )
#install(TARGETS oldmesher DESTINATION bin)

add_library(mesh3d
  ${SEMCOMMON}/h5helper.cpp
  ${SEMCOMMON}/read_unv.cpp
  ${SEMCOMMON}/read_input.c
  ${PRJ}/CSRC/sem_gll.cpp
  ${PRJ}/CSRC/mesh_h5_output.cpp
  ${PRJ}/CSRC/meshpart.cpp
  ${PRJ}/CSRC/mesh.cpp
  ${PRJ}/CSRC/earth_mesh.cpp
  ${PRJ}/CSRC/point3d.cpp
  ${PRJ}/CSRC/mesh_grid.cpp
  ${PRJ}/CSRC/mesh_common.cpp
  ${PRJ}/CSRC/reader_abaqus.cpp
  ${PRJ}/CSRC/reader_ideas.cpp
  ${SEM_SRC_C}
  )
#set_property(SOURCE ${SEMCOMMON}/read_unv.cpp ${PRJ}/CSRC/reader_ideas.cpp PROPERTY CXX_STANDARD 11) # set_property C++11 KO...
set_source_files_properties(${SEMCOMMON}/read_unv.cpp ${PRJ}/CSRC/reader_ideas.cpp PROPERTIES COMPILE_FLAGS "-std=c++11") # Bypass qui devrait marcher dans 99% des cas (gnu, intel)...

add_executable(sem_part_grid
  ${PRJ}/CSRC/main_grid.cpp
  )

add_executable(sem_part_h5
  ${PRJ}/CSRC/main_h5.cpp
  )

add_executable(mesher
  ${PRJ}/CSRC/main_mesher.cpp
  )
#set_property(SOURCE ${PRJ}/CSRC/main_mesher.cpp PROPERTY CXX_STANDARD 11) # set_property C++11 KO...
set_source_files_properties(${PRJ}/CSRC/main_mesher.cpp PROPERTIES COMPILE_FLAGS "-std=c++11") # Bypass qui devrait marcher dans 99% des cas (gnu, intel)...

include_directories(BEFORE ${PRJ}/CSRC ${PRJ}/metis/include)
include_directories(${BASE}/COMMON)
target_link_libraries(sem_part_grid mesh3d metis ${SEM_CPP_LIBS})
target_link_libraries(sem_part_h5 mesh3d metis ${SEM_CPP_LIBS})
target_link_libraries(mesher mesh3d metis ${SEM_CPP_LIBS})

install(TARGETS sem_part_grid sem_part_h5 mesher DESTINATION bin)
