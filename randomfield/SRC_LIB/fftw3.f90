module FFTW3
    use, intrinsic :: iso_c_binding
    !include 'fftw3.f'
    include 'fftw3-mpi.f03'
    !Attention there is a line "include 'fftw3-mpi.f03'" just before this one
end module FFTW3
