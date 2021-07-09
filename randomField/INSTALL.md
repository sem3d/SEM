# randomFields

## INSTALL

### External Libraries

`randomFields` is a MPI-based application to genrate 3D random fields. It exploits two main libraries: HDF5 and FFTW. In the following, Open MPI is considered as standard MPI implementation. 

Download [the latest OpenMPI release](https://www.open-mpi.org/software/ompi/v4.0) (tar.gz archive).

Download [the latest HDF5 release](https://www.hdfgroup.org/downloads/hdf5/source-code) (tar.gz archive). HDF5 must be installed by linking it to the available Open MPI distribution and specifying the installation prefix:

    ```
        tar -xvzf hdf5-version.tar.gz
        cd hdf5-version
        ./configure --prefix=$2 --enable-parallel --enable-fortran --enable-fortran2003 FC=mpif90 CC=mpicc
        make
        make check
        make install
    ```

Download [the latest FFTW3 release](http://fftw.org/download.html) (tar.gz archive). FFTW3 must be installed by linking it to the available Open MPI distribution and specifying the installation prefix:
    
    ```
        tar -xvzf fftw3-version.tar.gz
        cd fftw3-version
        ./configure --prefix=/path/to/install/fftw3 --enable-mpi --enable-openmp --enable-threads CC=gcc MPICC=mpicc F77=gfortran LDFLAGS=-L/path/to/openmpi/lib CPPFLAGS=-I/path/to/openmpi/include MPILIBS=-lmpi --disable-doc
        make
        make check
        make install        
    ```

### Install randomFields

Modify the headers of the template makefile provided `makefile_example` (Open MPI & mpif90 compiler) to correctly link the following libraries:

- MPI (Open MPI)

    ```
        LIBMPI = -L/path/to/local/mpi/lib -lmpi_cxx -lmpi -lmpi_mpifh -lmpi_usempif08 -lmpi_usempi_ignore_tkr
        INCLUDEMPI = -I/path/to/local/mpi/lib/include
    ```

 - HDF5

    ```
        LIBHDF5 = -L/path/to/local/hdf5/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -ldl -lm -lz 
        INCLUDEHDF5 = -I/path/to/local/hdf5/include
    ```
          
 - FFTW3

    ```
        LIBFFTW = -L/path/to/local/fftw3/lib -lfftw3_mpi -lfftw3 -lfftw3_threads
        INCLUDEFFTW = -I/path/to/local/fftw/include
    ```

The customized `makefile_example` must be moved into the build folder, located at the same level that `sem-ecp`:

```
    mkdir build_rf
    cp ./sem-ecp/randomFields/makefile_example ./build_rf/makefile
    cd build_rf
    make all
```