PROJECT(sem3d)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)

set(PRJ ${PROJECT_SOURCE_DIR})

set(SEM3D_SOURCES_F
    ${PRJ}/SRC/geometrical_prop.f90
    ${PRJ}/SRC/Modules/Element.f90
    ${PRJ}/SRC/Modules/Face.f90
    ${PRJ}/SRC/Modules/Vertex.f90
    ${PRJ}/SRC/Modules/Sources.f90
    ${PRJ}/SRC/Modules/TimeP.f90
    ${PRJ}/SRC/Modules/logical_input.f90
    ${PRJ}/SRC/Modules/Subdomains.f90
    ${PRJ}/SRC/Modules/Fault.f90
    ${PRJ}/SRC/Modules/Domain.f90
    ${PRJ}/SRC/Modules/Couplage.f90
    ${PRJ}/SRC/Modules/solid_fluid.f90
    ${PRJ}/SRC/orientation.f90
    ${PRJ}/SRC/assembly.f90
    ${PRJ}/SRC/Capteur.f90
    ${PRJ}/SRC/comm_utils.f90
    ${PRJ}/SRC/Newmark.f90
    ${PRJ}/SRC/solid_fluid_coupling.f90
    ${PRJ}/SRC/PML_def.f90
    ${PRJ}/SRC/SourcePosition.f90
    ${PRJ}/SRC/Modules/Edge.f90
    ${PRJ}/SRC/Modules/Comm.f90
    ${PRJ}/SRC/Modules/Neumann.f90
    ${PRJ}/SRC/Modules/PlaneW.f90
    ${PRJ}/SRC/Modules/Surface.f90
    ${PRJ}/SRC/Modules/tensor_util.f90
    ${PRJ}/SRC/allocate_domain.f90
    ${PRJ}/SRC/attenuation_aniso_update.f90
    ${PRJ}/SRC/attenuation_update.f90
    ${PRJ}/SRC/calcul_forces.f90
    ${PRJ}/SRC/calcul_forces_aniso.f90
    ${PRJ}/SRC/calcul_forces_aniso_att.f90
    ${PRJ}/SRC/calcul_forces_att.f90
    ${PRJ}/SRC/calcul_forces_fluid.f90
    ${PRJ}/SRC/compute_GLL.f90
    ${PRJ}/SRC/courant.f90
    ${PRJ}/SRC/deallocate_domain.f90
    ${PRJ}/SRC/define_arrays.f90
    ${PRJ}/SRC/define_Neumann_properties.f90
    ${PRJ}/SRC/define_planew_properties.f90
    ${PRJ}/SRC/double_couple.f90
    ${PRJ}/SRC/forces_int.f90
#    ${PRJ}/SRC/get_Displ_Edge2Elem.f90
#    ${PRJ}/SRC/get_Displ_Face2Elem.f90
#    ${PRJ}/SRC/get_Displ_Vertex2Elem.f90
#    ${PRJ}/SRC/get_Forces_Elem2Edge.f90
#    ${PRJ}/SRC/get_Forces_Elem2Face.f90
#    ${PRJ}/SRC/get_Forces_Elem2Vertex.f90
#    ${PRJ}/SRC/get_Mass_Elem2Edge.f90
#    ${PRJ}/SRC/get_Mass_Elem2Face.f90
#    ${PRJ}/SRC/get_Mass_Elem2Vertex.f90
#    ${PRJ}/SRC/getPMLpred_e2el.f90
#    ${PRJ}/SRC/getPMLpred_f2el.f90
#    ${PRJ}/SRC/getPMLpred_v2el.f90
    ${PRJ}/SRC/global_numbering.f90
    ${PRJ}/SRC/read_input.f90
    ${PRJ}/SRC/mesh3d.f90
    ${PRJ}/SRC/read_restart.f90
    ${PRJ}/SRC/save_checkpoint.f90
    ${PRJ}/SRC/second_fields.f90
    ${PRJ}/SRC/set_attenuation_param.f90
    ${PRJ}/SRC/shape8.f90
    ${PRJ}/SRC/shape27.f90
    ${PRJ}/SRC/drive_sem.f90
    ${PRJ}/SRC/ondelette.f90
    ${PRJ}/SRC/define_random.f90
    ${PRJ}/SRC/build_coord_mask.f90

#    ${PRJ}/SRC/get_VectProperty_Edge2Elem.f90
#    ${PRJ}/SRC/get_VectProperty_Elem2Edge.f90
#    ${PRJ}/SRC/get_VectProperty_Elem2Face.f90
#    ${PRJ}/SRC/get_VectProperty_Elem2Vertex.f90
#    ${PRJ}/SRC/get_VectProperty_Face2Elem.f90
#    ${PRJ}/SRC/get_VectProperty_Vertex2Elem.f90

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
    ${PRJ}/SRC/partial_deriv.f90

    ${PRJ}/SRC/snapshots.f90

	${PRJ}/SRC/earth_models/earth_transform.f90
	${PRJ}/SRC/earth_models/material_earthchunk.f90
	${PRJ}/SRC/earth_models/material_prem.f90
	${PRJ}/SRC/earth_models/model_earthchunk.f90
	${PRJ}/SRC/earth_models/model_prem.f90
	
)

#if(NOT USE_MKL)
#set(SEM3D_SOURCES_F
#    ${SEM3D_SOURCES_F}
#    ${PRJ}/SRC/BLAS/dgemm.f90
#    ${PRJ}/SRC/BLAS/lsame.f90
#    ${PRJ}/SRC/BLAS/xerbla.f90
#)
#endif(NOT USE_MKL)

set(SEM3D_SOURCES ${SEM3D_SOURCES_F}
    ${PRJ}/SRC/compute_constant_Q.c
    ${PRJ}/SRC/compute_iso_Q.c
    ${PRJ}/SRC/constant_Q2_sub.c
    ${PRJ}/SRC/iso_Q2_sub.c
    ${PRJ}/SRC/nrutil_alloc.c
    ${PRJ}/SRC/nrutil.c
)


if(OPT_MPI)
  #set(SEM_EXTRA_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
  # TODO: avec mpich MPI_LIBRARY contient libmpichcxx, qui a besoin des librairies de support
  # c++ et ne link pas semx.exe

  add_executable(sem3d_cpl ${PRJ}/SRC/main.F90 ${SEM3D_SOURCES})
  target_link_libraries(sem3d_cpl ${SEM_EXTRA_LIBS})
  set_target_properties(sem3d_cpl PROPERTIES 
                              COMPILE_DEFINITIONS COUPLAGE
                              COMPILE_FLAGS -I${PROJECT_BINARY_DIR}/cpl
                              COMPILE_FLAGS -I${SEMLIBDIR}
                              Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/cpl)
else(OPT_MPI)
  set(SEM3D_SOURCES ${SEM3D_SOURCES} ${FTM_SOURCES})
endif(OPT_MPI)

add_executable(sem3d.exe  ${PRJ}/SRC/main.F90 ${SEM3D_SOURCES})
target_link_libraries(sem3d.exe ${SEM_EXTRA_LIBS})
set_target_properties(sem3d.exe PROPERTIES 
                                COMPILE_FLAGS -I${PROJECT_BINARY_DIR}/ncpl
				COMPILE_FLAGS -I${SEMLIBDIR}
                                Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/ncpl)
install(TARGETS sem3d.exe RUNTIME DESTINATION bin)


add_executable(test_Q_coeffs ${PRJ}/TESTS/UNITS/test_Q_coeffs.c
      ${PRJ}/SRC/compute_constant_Q.c
      ${PRJ}/SRC/compute_iso_Q.c
      ${PRJ}/SRC/constant_Q2_sub.c
      ${PRJ}/SRC/iso_Q2_sub.c
      ${PRJ}/SRC/nrutil_alloc.c
      ${PRJ}/SRC/nrutil.c
)

target_link_libraries(test_Q_coeffs m)

add_executable(test_pointlist ${PRJ}/TESTS/UNITS/test_pointlist.f90
  ${SEMCOMMON}/pointlist.f90
)
