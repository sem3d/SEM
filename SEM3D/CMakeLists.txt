PROJECT(sem3d)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)

set(PRJ ${PROJECT_SOURCE_DIR})

set(SEM3D_SOURCES_F
    ${PRJ}/SRC/stat.f90
    ${PRJ}/SRC/geometrical_prop.f90
    ${PRJ}/SRC/Modules/Champs.f90
    ${PRJ}/SRC/Modules/Element.f90
    ${PRJ}/SRC/Modules/Face.f90
    ${PRJ}/SRC/Modules/Vertex.f90
    ${PRJ}/SRC/Modules/Sources.f90
    ${PRJ}/SRC/Modules/TimeP.f90
    ${PRJ}/SRC/Modules/logical_input.f90
    ${PRJ}/SRC/Modules/Subdomains.f90
    ${PRJ}/SRC/Modules/Fault.f90
    ${PRJ}/SRC/Modules/Interface.f90
    ${PRJ}/SRC/Modules/Domain.f90
    ${PRJ}/SRC/Modules/solid_fluid.f90
    ${PRJ}/SRC/orientation.f90
    ${PRJ}/SRC/Capteur.f90
    ${PRJ}/SRC/comm_utils.f90
    ${PRJ}/SRC/solid_fluid_coupling.f90
    ${PRJ}/SRC/SourcePosition.f90
    ${PRJ}/SRC/Modules/Edge.f90
    ${PRJ}/SRC/Modules/Comm.f90
    ${PRJ}/SRC/Modules/Neumann.f90
    ${PRJ}/SRC/Modules/PlaneW.f90
    ${PRJ}/SRC/Modules/Surface.f90
    ${PRJ}/SRC/Modules/tensor_util.f90
    ${PRJ}/SRC/Modules/Nonlinear.f90
    ${PRJ}/SRC/attenuation_aniso_update.f90
    ${PRJ}/SRC/attenuation_update.f90
    ${PRJ}/SRC/calcul_forces.f90
    ${PRJ}/SRC/calcul_forces_aniso.f90
    ${PRJ}/SRC/calcul_forces_aniso_att.f90
    ${PRJ}/SRC/calcul_forces_att.f90
    ${PRJ}/SRC/calcul_forces_fluid.f90
    ${PRJ}/SRC/compute_GLL.f90
    ${PRJ}/SRC/courant.f90
    ${PRJ}/SRC/define_arrays.f90
    ${PRJ}/SRC/define_Neumann_properties.f90
    ${PRJ}/SRC/define_planew_properties.f90
    ${PRJ}/SRC/double_couple.f90
    ${PRJ}/SRC/forces_int.f90
    ${PRJ}/SRC/global_numbering.f90
    ${PRJ}/SRC/mesh3d.f90
    ${PRJ}/SRC/read_input.f90
    ${PRJ}/SRC/read_restart.f90
    ${PRJ}/SRC/save_checkpoint.f90
    ${PRJ}/SRC/second_fields.f90
    ${PRJ}/SRC/set_attenuation_param.f90
    ${PRJ}/SRC/shape8.f90
    ${PRJ}/SRC/shape27.f90
    ${PRJ}/SRC/ondelette.f90
    ${PRJ}/SRC/define_random.f90
    ${PRJ}/SRC/build_prop_files.f90
    ${PRJ}/SRC/find_location.f90

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

set(SEM3D_SOURCES ${SEM3D_SOURCES_F})

# Sources qui dependent du flag COUPLAGE
set(SEM3D_SOURCES_CPL
    ${PRJ}/SRC/Modules/Couplage.f90
    ${PRJ}/SRC/Newmark.f90
    ${PRJ}/SRC/deallocate_domain.f90
    ${PRJ}/SRC/allocate_domain.f90
    ${PRJ}/SRC/drive_sem.f90
    ${PRJ}/SRC/main.F90
    )

add_library(sem3dcommon ${SEM3D_SOURCES})
set_target_properties(sem3dcommon PROPERTIES
				  COMPILE_FLAGS "-I${SEMLIBDIR}/sem3d -I${SEMLIBDIR}/common"
                                  Fortran_MODULE_DIRECTORY ${SEMLIBDIR}/sem3d)

if(OPT_MPI)
  #set(SEM_EXTRA_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
  # TODO: avec mpich MPI_LIBRARY contient libmpichcxx, qui a besoin des librairies de support
  # c++ et ne link pas semx.exe

#  add_executable(sem3d_cpl ${SEM3D_SOURCES_CPL})
#  target_link_libraries(sem3d_cpl sem3dcommon ${SEM_EXTRA_LIBS})
#  set_target_properties(sem3d_cpl PROPERTIES
#                                  COMPILE_DEFINITIONS COUPLAGE
#				  COMPILE_FLAGS "-I${SEMLIBDIR}/cpl3d -I${SEMLIBDIR}/sem3d -I${SEMLIBDIR}/common"
#				  Fortran_MODULE_DIRECTORY ${SEMLIBDIR}/cpl3d)
else(OPT_MPI)
  set(SEM3D_SOURCES ${SEM3D_SOURCES} ${FTM_SOURCES})
endif(OPT_MPI)

add_executable(sem3d.exe  ${SEM3D_SOURCES_CPL})
target_link_libraries(sem3d.exe sem3dcommon ${SEM_EXTRA_LIBS})
set_target_properties(sem3d.exe PROPERTIES 
				COMPILE_FLAGS "-I${SEMLIBDIR}/ncpl3d -I${SEMLIBDIR}/sem3d -I${SEMLIBDIR}/common"
                                Fortran_MODULE_DIRECTORY ${SEMLIBDIR}/ncpl3d)
install(TARGETS sem3d.exe RUNTIME DESTINATION bin)



add_executable(test_pointlist
  ${SEMCOMMON}/pointlist.f90
  ${PRJ}/TESTS/UNITS/test_pointlist.f90
)

add_executable(test_index
      ${PRJ}/SRC/indexation.f90
      ${PRJ}/TESTS/SRC/test_indexation.f90
      )