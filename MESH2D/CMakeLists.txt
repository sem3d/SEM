PROJECT(mesh2d)
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
ENABLE_LANGUAGE(Fortran CXX)
set(PRJ ${PROJECT_SOURCE_DIR})
set(BASE ${CMAKE_SOURCE_DIR})

set(MESHC_SRC
  ${BASE}/COMMON/h5helper.cpp
  ${SEMCOMMON}/read_unv.cpp
  ${PRJ}/SRC/main_mesh2d.cpp
  ${BASE}/COMMON/read_unv.cpp
)
#set_property(SOURCE ${SEMCOMMON}/read_unv.cpp ${PRJ}/SRC/main_mesh2d.cpp PROPERTY CXX_STANDARD 11) # set_property C++11 KO...
set_source_files_properties(${SEMCOMMON}/read_unv.cpp ${PRJ}/SRC/main_mesh2d.cpp PROPERTIES COMPILE_FLAGS "-std=c++11") # Bypass qui devrait marcher dans 99% des cas (gnu, intel)...

set(MESH_EXTRA_LIBS ${SEM_EXTRA_LIBS})

include_directories(BEFORE ${BASE}/MESH/metis/include)
include_directories(${BASE}/COMMON)
add_executable(mesh2dc ${MESHC_SRC})
target_link_libraries(mesh2dc metis ${MESH_EXTRA_LIBS})

install(TARGETS mesh2dc DESTINATION bin)
