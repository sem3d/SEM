PROJECT(sem3d)
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
ENABLE_LANGUAGE(Fortran)

set(PRJ ${PROJECT_SOURCE_DIR})

set(SEM3D_SOURCES_F
    ${PRJ}/SRC/stat.f90
    ${PRJ}/SRC/geometrical_prop.f90
    ${PRJ}/SRC/dombase.F90
    ${PRJ}/SRC/snapdata.F90
    ${PRJ}/SRC/Solid/champs_solid.f90
    ${PRJ}/SRC/SolidDG/champs_solid_dg.f90
    ${PRJ}/SRC/Fluid/champs_fluid.f90
    ${PRJ}/SRC/Modules/Element.f90
    ${PRJ}/SRC/Modules/Face.f90
    ${PRJ}/SRC/Modules/Vertex.f90
    ${PRJ}/SRC/Modules/Sources.f90
    ${PRJ}/SRC/Modules/Extended_source.f90
    ${PRJ}/SRC/Modules/TimeP.f90
    ${PRJ}/SRC/Modules/logical_input.f90
    ${PRJ}/SRC/Modules/Subdomains.f90
    ${PRJ}/SRC/Modules/Fault.f90
    ${PRJ}/SRC/Modules/Interface.f90
    ${PRJ}/SRC/Modules/Domain.f90
    ${PRJ}/SRC/Modules/solid_fluid.f90
    ${PRJ}/SRC/Modules/Nonlinear.f90
    ${PRJ}/SRC/indexation.f90
    ${PRJ}/SRC/global_numbering.f90
    ${PRJ}/SRC/partial_deriv.f90
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
#   ${PRJ}/SRC/Solid/attenuation_solid.F90
    ${PRJ}/SRC/Solid/calcul_forces.F90
    ${PRJ}/SRC/Solid/calcul_forces_nl.f90   
    ${PRJ}/SRC/Solid/calcul_forces_atn.F90
    # NOT WORKING
#    ${PRJ}/SRC/Solid/calcul_forces_atn_nl.F90
    ${PRJ}/SRC/Solid/dom_solid.F90
    ${PRJ}/SRC/Fluid/calcul_forces_fluid.F90
    ${PRJ}/SRC/Fluid/dom_fluid.F90
    ${PRJ}/SRC/SolidDG/calcul_forces_solid_dg.F90
    ${PRJ}/SRC/SolidDG/dom_solid_dg.F90
    ${PRJ}/SRC/compute_GLL.f90
    ${PRJ}/SRC/courant.f90
    ${PRJ}/SRC/define_arrays.F90
    ${PRJ}/SRC/define_Neumann_properties.f90
    ${PRJ}/SRC/double_couple.f90
    ${PRJ}/SRC/mesh3d.f90
    ${PRJ}/SRC/read_input.f90
    ${PRJ}/SRC/read_restart.f90
    ${PRJ}/SRC/save_checkpoint.f90
    ${PRJ}/SRC/set_attenuation_param.F90
    ${PRJ}/SRC/shape8.F90
    ${PRJ}/SRC/shape27.F90
    ${PRJ}/SRC/ondelette.f90
    ${PRJ}/SRC/find_location.f90
    ${PRJ}/SRC/charFunctions.f90

    ${PRJ}/SRC/algo_comm.f90
    ${PRJ}/SRC/source_excit.F90
#    ${PRJ}/SRC/boundary_conditions_init.f90

    ${PRJ}/SRC/snapshots.F90

    ${PRJ}/SRC/build_prop_files.F90

    ${PRJ}/SRC/earth_models/earth_transform.f90
    ${PRJ}/SRC/earth_models/material_earthchunk.F90
    ${PRJ}/SRC/earth_models/material_prem.F90
    ${PRJ}/SRC/earth_models/model_earthchunk.f90
    ${PRJ}/SRC/earth_models/model_prem.f90

    ${PRJ}/SRC/SourceMT/parameters.f90
    ${PRJ}/SRC/SourceMT/mathfval.f90
    ${PRJ}/SRC/SourceMT/fvalue.f90
    ${PRJ}/SRC/SourceMT/func3D.f90
    ${PRJ}/SRC/SourceMT/func2D.f90
    ${PRJ}/SRC/SourceMT/func1D.f90
    ${PRJ}/SRC/SourceMT/fparser.f90
    ${PRJ}/SRC/SourceMT/ErrorWarning.f90
    ${PRJ}/SRC/SourceMT/Dimension.f90
    ${PRJ}/SRC/SourceMT/create_sem_surface_input.f90
    ${PRJ}/SRC/SourceMT/compute_surface_BC.f90
    ${PRJ}/SRC/SourceMT/DimensionalShape.f90
    ${PRJ}/SRC/SourceMT/selectsurfload.f90 
    ${PRJ}/SRC/Modules/Bspline.f90
    ${PRJ}/SRC/Modules/Mirror.f90

)

if (OPT_CPML)
  set(SEM3D_SOURCES_F ${SEM3D_SOURCES_F}
      ${PRJ}/SRC/FluidCPML/champs_fluidcpml.f90
      ${PRJ}/SRC/SolidCPML/champs_solidcpml.f90
      ${PRJ}/SRC/pml.F90
      ${PRJ}/SRC/FluidCPML/calcul_forces_fluidpml.F90
      ${PRJ}/SRC/FluidCPML/dom_fluidcpml.F90
      ${PRJ}/SRC/SolidCPML/calcul_forces_solidpml.F90
      ${PRJ}/SRC/SolidCPML/dom_solidcpml.F90)
else(OPT_CPML)
  set(SEM3D_SOURCES_F ${SEM3D_SOURCES_F}
      ${PRJ}/SRC/FluidPML/champs_fluidpml.f90
      ${PRJ}/SRC/SolidPML/champs_solidpml.f90
      ${PRJ}/SRC/pml.F90
      ${PRJ}/SRC/FluidPML/dom_fluidpml.F90
      ${PRJ}/SRC/SolidPML/dom_solidpml.F90)
endif(OPT_CPML)

#if(NOT USE_MKL)
#set(SEM3D_SOURCES_F
#    ${SEM3D_SOURCES_F}
#    ${PRJ}/SRC/BLAS/dgemm.f90
#    ${PRJ}/SRC/BLAS/lsame.f90
#    ${PRJ}/SRC/BLAS/xerbla.f90
#)
#endif(NOT USE_MKL)

set(SEM3D_SOURCES ${SEM3D_SOURCES_F})

set(SEM3D_SOURCES_MAIN
    ${PRJ}/SRC/Newmark.f90
    ${PRJ}/SRC/deallocate_domain.f90
    ${PRJ}/SRC/allocate_domain.f90
    ${PRJ}/SRC/drive_sem.f90
    ${PRJ}/SRC/main.F90
    )

if(OPT_MPI)
  #set(SEM_EXTRA_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
  # TODO: avec mpich MPI_LIBRARY contient libmpichcxx, qui a besoin des librairies de support
  # c++ et ne link pas semx.exe
else(OPT_MPI)
  set(SEM3D_SOURCES ${SEM3D_SOURCES} ${FTM_SOURCES})
endif(OPT_MPI)

include_directories(${PRJ}/SRC)

add_executable(sem3d.exe  ${SEM3D_SOURCES} ${SEM3D_SOURCES_MAIN} )
target_link_libraries(sem3d.exe ${SEM_EXTRA_LIBS})
set_target_properties(sem3d.exe PROPERTIES
  COMPILE_FLAGS "-I${SEMLIBDIR}/sem3d -I${SEMLIBDIR}/common"
  Fortran_MODULE_DIRECTORY ${SEMLIBDIR}/sem3d)
install(TARGETS sem3d.exe RUNTIME DESTINATION bin)


add_executable(test_pointlist
  ${SEMCOMMON}/pointlist.f90
  ${PRJ}/TESTS/UNITS/test_pointlist.f90
)

add_executable(test_index
      ${PRJ}/SRC/indexation.f90
      ${PRJ}/TESTS/SRC/test_indexation.f90
      )
