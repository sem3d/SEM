cmake_minimum_required(VERSION 3.10)
PROJECT(randomField Fortran C)

set(RFDIR ${PROJECT_SOURCE_DIR})

file(GLOB DCDFLIB_SRC "${RFDIR}/dcdflib/src/*.f")

# Random Field
set(RF_SRC
  ${RFDIR}/SRC_EXE/main_RandomField.f90
  ${RFDIR}/SRC_LIB/Rosemblat.f90
  ${RFDIR}/SRC_LIB/randomFieldND.f90
  ${RFDIR}/SRC_LIB/spectra_RF.f90
  ${RFDIR}/SRC_LIB/writeResultFile_RF.f90
  ${RFDIR}/SRC_LIB/mesh_RF.f90
  ${RFDIR}/SRC_LIB/type_MESH.f90
  ${RFDIR}/SRC_LIB/readUNV_RF.f90
  ${RFDIR}/SRC_LIB/type_RF.f90
  ${RFDIR}/SRC_LIB/topography_RF.f90
  ${RFDIR}/SRC_LIB/localization_RF.f90
  ${RFDIR}/SRC_LIB/type_inputRF.f90
  ${RFDIR}/SRC_LIB/sample_RF.f90
  ${RFDIR}/SRC_LIB/calls_RF.f90
  ${DCDFLIB_SRC}
  ${RFDIR}/RANDLIB/ranlib.f90
  ${RFDIR}/RANDLIB/rnglib.f90
  )

set(STAT_SRC
  ${RFDIR}/SRC_LIB/statistics_RF.f90
  ${RFDIR}/SRC_LIB/type_STAT.f90
  ${RFDIR}/SRC_EXE/main_Stat.f90
  )
set(COMM_SRC
  ${RFDIR}/SRC_LIB/special_functions.f90
  ${RFDIR}/SRC_LIB/hdf5_RF.f90
  ${RFDIR}/SRC_LIB/readFile_RF.f90
  ${RFDIR}/SRC_LIB/systemUt_RF.f90
  ${RFDIR}/SRC_LIB/charFunctions.f90
  ${RFDIR}/SRC_LIB/constants_RF.f90
  ${RFDIR}/SRC_LIB/common_variables_RF.f90
  ${RFDIR}/SRC_LIB/displayCarvalhol.f90
  ${RFDIR}/SRC_LIB/write_Log_File.f90
  ${RFDIR}/SRC_LIB/math_RF.f90
  ${RFDIR}/SRC_LIB/fftw3.f90
  ${RFDIR}/BLAS/dgemm.f90
  ${RFDIR}/BLAS/lsame.f90
  ${RFDIR}/BLAS/xerbla.f90
  ${RFDIR}/BLAS/dger.f90
  ${RFDIR}/SRC_LIB/hdf5_compat.c
  )

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
link_directories(/opt/homebrew/lib)

find_package(FFTW3 REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(/opt/homebrew/include)

message("FFTW3: ${FFTW3_LIBRARIES}")
include_directories(${FFTW3_INCLUDE_DIR})
add_executable(randomField ${RF_SRC} ${COMM_SRC})
add_executable(statistics ${STAT_SRC} ${COMM_SRC})
target_link_libraries(randomField ${FFTW3_LIBRARIES} fftw3_mpi hdf5_hl_fortran hdf5_hl hdf5_fortran hdf5 ${HDF5_LIBRARIES} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
target_link_libraries(statistics ${FFTW3_LIBRARIES} fftw3_mpi hdf5_hl_fortran hdf5_hl hdf5_fortran hdf5 ${HDF5_LIBRARIES} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
