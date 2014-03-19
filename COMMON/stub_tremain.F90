


subroutine tremain( remtime )
    interface
       subroutine tremain_c(remtime) bind(c)
           use iso_c_binding
           real(C_DOUBLE), intent(out) :: remtime
       end subroutine tremain_c
    end interface
    real(kind=8), intent(out) :: remtime

    remtime = 1e12
    call tremain_c(remtime)

end subroutine tremain
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
