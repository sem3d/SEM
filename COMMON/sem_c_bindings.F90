!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module sem_c_bindings
    use iso_c_binding

    public :: sem_mkdir, strlen, sem_check_file

    interface
       function strlen(s) bind(c)
           use iso_c_binding
           type(c_ptr), intent(in), value :: s
           integer(C_INT) :: strlen
       end function strlen

       function sem_mkdir_c(dirname) bind(c)
           use iso_c_binding
           character(C_CHAR), dimension(*) :: dirname
           integer(C_INT) :: sem_mkdir_c
       end function sem_mkdir_c

       function sem_check_file_c(filename) bind(c)
           use iso_c_binding
           character(C_CHAR), dimension(*) :: filename
           integer(C_INT) :: sem_check_file_c
       end function sem_check_file_c

    end interface

contains

    function sem_mkdir(dirname)
        use iso_c_binding
        character(kind=C_CHAR,len=*) :: dirname
        integer(C_INT) :: sem_mkdir
        character(C_CHAR), dimension(:), allocatable :: c_name
        integer :: i, l

        l = len_trim(dirname)
        allocate(c_name(1:l+1))
        do i=1,l
            c_name(i) = dirname(i:i)
        enddo
        c_name(l+1) = C_NULL_CHAR
        sem_mkdir = sem_mkdir_c(c_name)
        deallocate(c_name)
    end function sem_mkdir

    function sem_check_file(filename)
        use iso_c_binding
        character(kind=C_CHAR,len=*) :: filename
        integer(C_INT) :: sem_check_file
        character(C_CHAR), dimension(:), allocatable :: c_name
        integer :: i, l

        l = len_trim(filename)
        allocate(c_name(1:l+1))
        do i=1,l
            c_name(i) = filename(i:i)
        enddo
        c_name(l+1) = C_NULL_CHAR
        sem_check_file = sem_check_file_c(c_name)
        deallocate(c_name)
    end function sem_check_file


end module sem_c_bindings

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
