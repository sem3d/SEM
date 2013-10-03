PROJECT(sem2d)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_LANGUAGE(Fortran)
set(PRJ ${PROJECT_SOURCE_DIR})

set(SEM2D_SOURCES_F
    ${SEM_SRC}
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
    ${PRJ}/SRC/Capteur.F90
    ${PRJ}/SRC/Couplage.F90
    ${PRJ}/SRC/Newmark.F90
    ${PRJ}/SRC/PML_def.F90
    ${PRJ}/SRC/ajout.F90
    ${PRJ}/SRC/allocate_domain.F90
    ${PRJ}/SRC/compute_GLL.F90
    ${PRJ}/SRC/convertirUnv.F90
    ${PRJ}/SRC/courant.F90
    ${PRJ}/SRC/define_arr.F90
    ${PRJ}/SRC/define_fault_properties.F90
    ${PRJ}/SRC/dumptrace.F90
    ${PRJ}/SRC/getDispl_fvel.F90
    ${PRJ}/SRC/getPMLpred_fv2el.F90
    ${PRJ}/SRC/getdumpmass_e2f.F90
    ${PRJ}/SRC/getif_el2f.F90
    ${PRJ}/SRC/getiv_e2f.F90
    ${PRJ}/SRC/getmass_e2f.F90
    ${PRJ}/SRC/global_numbering.F90
#    ${PRJ}/SRC/lister_faces_sem.f90
#    ${PRJ}/SRC/plot_grid.F90
#    ${PRJ}/SRC/postraitement.F90
    ${PRJ}/SRC/read_input.F90
    ${PRJ}/SRC/read_mesh.F90
    ${PRJ}/SRC/read_restart.F90
    ${PRJ}/SRC/receiver_pos.F90
    ${PRJ}/SRC/save_checkpoint.F90
    ${PRJ}/SRC/save_deformation.F90
    ${PRJ}/SRC/save_fault_trace.F90
    ${PRJ}/SRC/save_trace.F90
    ${PRJ}/SRC/save_vorticity.F90
    ${PRJ}/SRC/savefield.F90
    ${PRJ}/SRC/sem.F90
    ${PRJ}/SRC/shape4.F90
    ${PRJ}/SRC/shape8.F90
    ${PRJ}/SRC/smoothing_exp.F90
    ${PRJ}/SRC/source_pos.F90
    ${PRJ}/SRC/unv.F90
    ${PRJ}/SRC/verify_in_quad.F90
    ${PRJ}/SRC/wall_transfer.F90
    ${PRJ}/SRC/main.F90
)

set(SEM2D_SOURCES ${SEM2D_SOURCES_F}
  ${SEM_SRC_C}
)


if (OPT_COUPLAGE)
add_library(sem2d STATIC ${SEM2D_SOURCES})
set_target_properties(sem2d PROPERTIES COMPILE_DEFINITIONS COUPLAGE)
install(TARGETS sem2d ARCHIVE DESTINATION lib)
endif (OPT_COUPLAGE)

# NB: on peut compiler sem2d.exe avec OPT_COUPLAGE mais on va compiler
# les sources 2 fois (une fois avec -DCOUPLAGE, une fois sans)
if(NOT OPT_COUPLAGE)
add_executable(sem2d.exe  ${SEM2D_SOURCES})
target_link_libraries(sem2d.exe ${CEAUSER} ${SEM_EXTRA_LIBS})
install(TARGETS sem2d.exe DESTINATION bin)
endif(NOT OPT_COUPLAGE)

