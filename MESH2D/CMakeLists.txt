PROJECT(mesh2d)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran CXX)
set(PRJ ${PROJECT_SOURCE_DIR})
set(BASE ${CMAKE_SOURCE_DIR})


set(MESHC_SRC
  ${BASE}/COMMON/h5helper.cpp
  ${PRJ}/SRC/main_mesh2d.cpp
)

set(MESH_EXTRA_LIBS ${SEM_EXTRA_LIBS})

include_directories(BEFORE ${BASE}/MESH/metis/include)
include_directories(${BASE}/COMMON)
add_executable(mesh2dc ${MESHC_SRC})
target_link_libraries(mesh2dc metis ${MESH_EXTRA_LIBS})

install(TARGETS mesh2dc DESTINATION bin)
