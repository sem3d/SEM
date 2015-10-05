PROJECT(sem2d)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)

set(PRJ ${PROJECT_SOURCE_DIR})

set(SEM2D_SOURCES
    ${PRJ}/SRC/Element.F90
    ${PRJ}/SRC/Face.F90
    ${PRJ}/SRC/Vertex.F90
    ${PRJ}/SRC/Sources.F90
    ${PRJ}/SRC/TimeP.F90
    ${PRJ}/SRC/logical_input.F90
    ${PRJ}/SRC/Subdomains.F90
    ${PRJ}/SRC/Receivers.F90
    ${PRJ}/SRC/Fault.F90
    ${PRJ}/SRC/mpi_list.F90
    ${PRJ}/SRC/communicating_wall.F90
    ${PRJ}/SRC/Domain.F90
    ${PRJ}/SRC/orientation.F90
    ${PRJ}/SRC/Capteur.F90
    ${PRJ}/SRC/Couplage.F90
    ${PRJ}/SRC/Newmark.F90
    ${PRJ}/SRC/PML_def.F90
    ${PRJ}/SRC/allocate_domain.F90
    ${PRJ}/SRC/compute_GLL.F90
    ${PRJ}/SRC/courant.F90
    ${PRJ}/SRC/define_arr.F90
    ${PRJ}/SRC/define_fault_properties.F90
    ${PRJ}/SRC/getDispl_fvel.F90
    ${PRJ}/SRC/getPMLpred_fv2el.F90
    ${PRJ}/SRC/getdumpmass_e2f.F90
    ${PRJ}/SRC/getif_el2f.F90
    ${PRJ}/SRC/getiv_e2f.F90
    ${PRJ}/SRC/getmass_e2f.F90
    ${PRJ}/SRC/global_energy.F90
    ${PRJ}/SRC/global_numbering.F90
    ${PRJ}/SRC/snapshots.F90
    ${PRJ}/SRC/read_input.F90
    ${PRJ}/SRC/shape4.F90
    ${PRJ}/SRC/shape8.F90
    ${PRJ}/SRC/find_location.F90
    ${PRJ}/SRC/receiver_pos.F90
    ${PRJ}/SRC/read_mesh.F90
    ${PRJ}/SRC/read_mesh_h5.F90
    ${PRJ}/SRC/read_restart.F90
    ${PRJ}/SRC/save_checkpoint.F90
    ${PRJ}/SRC/save_deformation.F90
    ${PRJ}/SRC/save_fault_trace.F90
    ${PRJ}/SRC/save_vorticity.F90
    ${PRJ}/SRC/savefield.F90
    ${PRJ}/SRC/smoothing_exp.F90
    ${PRJ}/SRC/source_pos.F90
    ${PRJ}/SRC/wall_transfer.F90
)

add_library(sem2dcommon ${SEM2D_SOURCES})
set_target_properties(sem2dcommon PROPERTIES
				  COMPILE_FLAGS "-I${SEMLIBDIR}/common -I${SEMLIBDIR}/sem2d"
                                  Fortran_MODULE_DIRECTORY ${SEMLIBDIR}/sem2d)


if(OPT_MPI)
  add_executable(sem2d_cpl ${PRJ}/SRC/main.F90 ${SEM2D_SOURCES})
  target_link_libraries(sem2d_cpl sem2dcommon ${SEM_EXTRA_LIBS})
  set_target_properties(sem2d_cpl PROPERTIES
                                  COMPILE_DEFINITIONS COUPLAGE
				  COMPILE_FLAGS "-I${SEMLIBDIR}/cpl2d -I${SEMLIBDIR}/common -I${SEMLIBDIR}/sem2d"
                              	  Fortran_MODULE_DIRECTORY ${SEMLIBDIR}/cpl2d)
  install(TARGETS sem2d_cpl DESTINATION bin)
else(OPT_MPI)
  set(SEM2D_SOURCES ${SEM2D_SOURCES} ${FTM_SOURCES})
endif(OPT_MPI)

add_executable(sem2d.exe ${PRJ}/SRC/main.F90 ${SEM2D_SOURCES})
target_link_libraries(sem2d.exe sem2dcommon ${SEM_EXTRA_LIBS})

set_target_properties(sem2d.exe PROPERTIES
				COMPILE_FLAGS "-I${SEMLIBDIR}/ncpl2d -I${SEMLIBDIR}/common -I${SEMLIBDIR}/sem2d"
                              	Fortran_MODULE_DIRECTORY ${SEMLIBDIR}/ncpl2d)

install(TARGETS sem2d.exe DESTINATION bin)

