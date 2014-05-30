PROJECT(mesh2d)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ENABLE_LANGUAGE(Fortran CXX)
set(PRJ ${PROJECT_SOURCE_DIR})
set(BASE ${CMAKE_SOURCE_DIR})


set(MESHC_SRC
  ${PRJ}/SRC/main_mesh2d.cpp
)

set(MESH_EXTRA_LIBS "")

if(OPT_MPI OR HDF5_IS_PARALLEL)
set(SEM_EXTRA_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
endif(OPT_MPI OR HDF5_IS_PARALLEL)

set(MESH_EXTRA_LIBS ${MESH_EXTRA_LIBS} ${HDF5_LIBRARIES})

include_directories(${BASE}/MESH/metis/include)

add_executable(mesh2dc ${MESHC_SRC})
target_link_libraries(mesh2dc metis ${MESH_EXTRA_LIBS})

install(TARGETS mesh2dc DESTINATION bin)
