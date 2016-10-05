module constants_RF

    implicit none

    double precision, parameter :: PI = 3.1415926535898d0;
    double precision, parameter :: SQRT_2PI = 2.50662827463d0;
    double precision, parameter :: ZERO = 0d0;
    double precision, parameter :: TOLERANCE = 0.000000001
    double precision, parameter :: MIN_DOUBLE = -2.0D+307
    double precision, parameter :: MAX_DOUBLE = 2.0D+307
	!METHOD
    integer, parameter :: ISOTROPIC = 1, &
                          SHINOZUKA = 2, &
                          RANDOMIZATION = 3, &
                          FFT = 4
    !Correlation Model
    integer, parameter :: cm_GAUSSIAN = 1, &
                          cm_COS = 2
    !First-order Marginal Density
    integer, parameter :: fom_GAUSSIAN = 1, &
                          fom_LOGNORMAL = 2
    !Mesh Mode
    integer, parameter :: msh_AUTO = 1, msh_UNV = 2

    integer, parameter :: SCREEN=6
    integer, parameter :: buf_RF=1024 !Buffer for text

    !Output
    logical :: DIRECT_OUT=.true. !When true output matrix is XYZ, when false the output is ZYX

    !APPLICATION
    integer, parameter :: NATIVE = 1
    integer, parameter :: SEM = 2

    integer :: TESTRANK = 0 !ONLY FOR TESTS
    character (len=30), parameter :: mesh_input = "mesh_input"
    character (len=30), parameter :: gen_input  = "gen_input"
    character (len=30), parameter :: test_input = "test_input"
    !character (len=30), parameter :: unv_input = "unv_files/Maroua1.unv "!"unv_files/Luciano_Cube.unv"

end module constants_RF
