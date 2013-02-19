PROJECT(sem3d)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)
set(PRJ ${PROJECT_SOURCE_DIR})


set(SEM3D_SOURCES_F
    ${SEM_SRC}
    ${PRJ}/SRC/geometrical_prop.f90
    ${PRJ}/SRC/Modules/Element.f90
    ${PRJ}/SRC/Modules/Face.f90
    ${PRJ}/SRC/Modules/Vertex.f90
    ${PRJ}/SRC/Modules/Sources.f90
    ${PRJ}/SRC/Modules/TimeP.f90
    ${PRJ}/SRC/Modules/logical_input.f90
    ${PRJ}/SRC/Modules/Subdomains.f90
    ${PRJ}/SRC/Modules/Receivers.f90
    ${PRJ}/SRC/Modules/Fault.f90
    ${PRJ}/SRC/Modules/Domain.f90
    ${PRJ}/SRC/Modules/Couplage.f90
    ${PRJ}/SRC/Modules/solid_fluid.f90
    ${PRJ}/SRC/Capteur.f90
    ${PRJ}/SRC/comm_utils.f90
    ${PRJ}/SRC/Newmark.f90
    ${PRJ}/SRC/solid_fluid_coupling.f90
    ${PRJ}/SRC/PML_def.f90
    ${PRJ}/SRC/ReceiverPosition.f90
    ${PRJ}/SRC/SourcePosition.f90
    ${PRJ}/SRC/Modules/Edge.f90
    ${PRJ}/SRC/Modules/Comm.f90
    ${PRJ}/SRC/Modules/Neumann.f90
    ${PRJ}/SRC/Modules/PlaneW.f90
    ${PRJ}/SRC/Modules/Surface.f90
    ${PRJ}/SRC/allocate_domain.f90
    ${PRJ}/SRC/attenuation_aniso_update.f90
    ${PRJ}/SRC/attenuation_update.f90
    ${PRJ}/SRC/calcul_forces.f90
    ${PRJ}/SRC/calcul_forces_aniso.f90
    ${PRJ}/SRC/calcul_forces_aniso_att.f90
    ${PRJ}/SRC/calcul_forces_att.f90
    ${PRJ}/SRC/compute_GLL.f90
    ${PRJ}/SRC/courant.f90
    ${PRJ}/SRC/deallocate_domain.f90
    ${PRJ}/SRC/define_arrays.f90
    ${PRJ}/SRC/define_Neumann_properties.f90
    ${PRJ}/SRC/define_planew_properties.f90
    ${PRJ}/SRC/double_couple.f90
    ${PRJ}/SRC/forces_int.f90
    ${PRJ}/SRC/get_Displ_Edge2Elem.f90
    ${PRJ}/SRC/get_Displ_Face2Elem.f90
    ${PRJ}/SRC/get_Displ_Vertex2Elem.f90
    ${PRJ}/SRC/get_Forces_Elem2Edge.f90
    ${PRJ}/SRC/get_Forces_Elem2Face.f90
    ${PRJ}/SRC/get_Forces_Elem2Vertex.f90
    ${PRJ}/SRC/get_Mass_Elem2Edge.f90
    ${PRJ}/SRC/get_Mass_Elem2Face.f90
    ${PRJ}/SRC/get_Mass_Elem2Vertex.f90
    ${PRJ}/SRC/getPMLpred_e2el.f90
    ${PRJ}/SRC/getPMLpred_f2el.f90
    ${PRJ}/SRC/getPMLpred_v2el.f90
    ${PRJ}/SRC/global_numbering.f90
    ${PRJ}/SRC/read_input.f90
    ${PRJ}/SRC/mesh3d.f90
    ${PRJ}/SRC/read_restart.f90
    ${PRJ}/SRC/save_checkpoint.f90
    ${PRJ}/SRC/savetrace.f90
    ${PRJ}/SRC/set_attenuation_aniso_param.f90
    ${PRJ}/SRC/set_attenuation_param.f90
    ${PRJ}/SRC/shape8.f90
    ${PRJ}/SRC/shape27.f90
    ${PRJ}/SRC/drive_sem.f90
    ${PRJ}/SRC/BLAS/dgemm.f90
    ${PRJ}/SRC/BLAS/lsame.f90
    ${PRJ}/SRC/BLAS/xerbla.f90
    ${PRJ}/SRC/ondelette.f90

    ${PRJ}/SRC/get_VectProperty_Edge2Elem.f90
    ${PRJ}/SRC/get_VectProperty_Elem2Edge.f90
    ${PRJ}/SRC/get_VectProperty_Elem2Face.f90
    ${PRJ}/SRC/get_VectProperty_Elem2Vertex.f90
    ${PRJ}/SRC/get_VectProperty_Face2Elem.f90
    ${PRJ}/SRC/get_VectProperty_Vertex2Elem.f90

    ${PRJ}/SRC/get_ScalarProperty_Edge2Elem.f90
    ${PRJ}/SRC/get_ScalarProperty_Elem2Edge.f90
    ${PRJ}/SRC/get_ScalarProperty_Elem2Face.f90
    ${PRJ}/SRC/get_ScalarProperty_Elem2Vertex.f90
    ${PRJ}/SRC/get_ScalarProperty_Face2Elem.f90
    ${PRJ}/SRC/get_ScalarProperty_Vertex2Elem.f90

    ${PRJ}/SRC/indexation.f90
    ${PRJ}/SRC/algo_comm.f90
    ${PRJ}/SRC/source_excit.f90
    ${PRJ}/SRC/boundary_conditions_init.f90

    ${PRJ}/SRC/snapshots.f90
)

set(SEM3D_SOURCES ${SEM3D_SOURCES_F}
    ${SEM_SRC_C}
    ${PRJ}/SRC/compute_constant_Q.c
    ${PRJ}/SRC/compute_iso_Q.c
    ${PRJ}/SRC/constant_Q2_sub.c
    ${PRJ}/SRC/iso_Q2_sub.c
    ${PRJ}/SRC/nrutil.c
)


# On veut compiler les .f90 avec le preprocesseur cpp
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set_source_files_properties(${SEM3D_SOURCES_F} ${FTM_SOURCES} PROPERTIES COMPILE_FLAGS "-fpp")
else (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set_source_files_properties(${SEM3D_SOURCES_F} ${FTM_SOURCES} PROPERTIES COMPILE_FLAGS "-x f95-cpp-input")
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")


set(SEM_EXTRA_LIBS "")
if(OPT_MPI)
#set(SEM_EXTRA_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
# TODO: avec mpich MPI_LIBRARY contient libmpichcxx, qui a besoin des librairies de support
# c++ et ne link pas semx.exe
set(SEM_EXTRA_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
else(OPT_MPI)
set(SEM3D_SOURCES ${SEM3D_SOURCES}
  ${FTM_SOURCES})
endif(OPT_MPI)

if (OPT_HDF5)
set(SEM_EXTRA_LIBS ${SEM_EXTRA_LIBS} hdf5_fortran)
set(SEM3D_SOURCES ${SEM3D_SOURCES} ${SEMCOMMON}/sem_hdf5.f90)
endif(OPT_HDF5)

if (OPT_COUPLAGE)
add_library(sem3d STATIC ${SEM3D_SOURCES})
set_target_properties(sem3d PROPERTIES COMPILE_DEFINITIONS COUPLAGE)
install(TARGETS sem3d ARCHIVE DESTINATION lib)
endif (OPT_COUPLAGE)

# NB: on peut compiler sem3d.exe avec OPT_COUPLAGE mais on va compiler
# les sources 2 fois (une fois avec -DCOUPLAGE, une fois sans)
if(NOT OPT_COUPLAGE)
add_executable(sem3d.exe  ${PRJ}/SRC/main.F90 ${SEM3D_SOURCES})
target_link_libraries(sem3d.exe ${CEAUSER} ${SEM_EXTRA_LIBS})
install(TARGETS sem3d.exe RUNTIME DESTINATION bin)
endif(NOT OPT_COUPLAGE)


add_executable(test_Q_coeffs ${PRJ}/TESTS/UNITS/test_Q_coeffs.c
      ${PRJ}/SRC/compute_constant_Q.c
    ${PRJ}/SRC/compute_iso_Q.c
    ${PRJ}/SRC/constant_Q2_sub.c
    ${PRJ}/SRC/iso_Q2_sub.c
    ${PRJ}/SRC/nrutil.c
)

target_link_libraries(test_Q_coeffs m)
