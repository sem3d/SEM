PROJECT(mesh2d)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran CXX)
set(PRJ ${PROJECT_SOURCE_DIR})
set(BASE ${CMAKE_SOURCE_DIR})


set(MESHC_SRC
  ${PRJ}/SRC/main_mesh2d.cpp
)

set(MESH_EXTRA_LIBS "")

if(OPT_MPI)
set(SEM_EXTRA_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
else(OPT_MPI)
endif(OPT_MPI)

add_definitions(-DUSE_HDF5)
include_directories(${HDF5_INCLUDE_DIR_FORTRAN})
#set(MESH_EXTRA_LIBS ${MESH_EXTRA_LIBS} hdf5_fortran)


add_executable(mesh2dc ${MESHC_SRC})
target_link_libraries(mesh2dc ${MESH_EXTRA_LIBS})
include_directories(${BASE}/MESH/metis/include ${HDF5_INCLUDE_DIR})
target_link_libraries(mesh2dc metis ${HDF5_LIBRARIES})

install(TARGETS mesh2dc DESTINATION bin)
