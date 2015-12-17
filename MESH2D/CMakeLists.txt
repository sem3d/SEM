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

# A reprendre plus tard : il faut imposer C++11 pour utiliser read_unv.cpp
set_property(TARGET mesh2dc PROPERTY CXX_STANDARD 11)          # read_unv.cpp needs C++11 (tuple) - Ne marche pas ?!
set_property(TARGET mesh2dc PROPERTY CXX_STANDARD_REQUIRED ON) # read_unv.cpp needs C++11 (tuple) - Ne marche pas ?!

install(TARGETS mesh2dc DESTINATION bin)
