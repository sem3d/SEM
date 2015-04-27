!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file pol_force.F90
!!\brief Contient la subroutine pol_force.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Assure l'évaluation de la valeur recherchée en x pour la fonction au point de gauss k.
!!
!! \param integer, intent (IN) n
!! \param integer, intent (IN) k
!! \param real, dimension (0:n-1), intent (IN) GLLc
!! \param real,intent (IN) x
!! \param real, intent (OUT) y
!<


subroutine pol_force (n,GLLc,k,x,y)
    !  n nombre de point de gauss
    !  GLLc valeur des points de gauss
    !  k indice de la fonction
    !  x  abscisse
    !  y  valeur recherche en x pour la fonction au point de gauss k
    implicit none
    integer, intent (IN) :: n,k
    real, dimension  (0:n-1), intent (IN) :: GLLc
    real,intent (IN) :: x
    real, intent (OUT) :: y

    y = 1
    if (n ==0 ) then
        write (*,*) "Bad number n = 0. It should be a constant!"
        stop
    endif
    if (n ==1 ) return
    ! ################################################################
    !  recherche de la valeur
    y = 0.
    if ( k .eq. 0 ) then
        if ( x .lt. GLLc(1) ) then
            y = (x -GLLc(1))/(GLLc(0) - GLLc(1))
        endif
    else
        if ( k .eq. n-1 ) then
            if ( x .gt. GLLc(n-2) ) then
                y = (x -GLLc(n-2))/(GLLc(n-1) - GLLc(n-2))
            endif
        else
            if ( x .gt. GLLc(k-1) .and. x .lt. Gllc(k+1) ) then
                if (  x .lt. GLLc(k) ) then
                    y = (x - GLLc(k-1))/(GLLc(k)-GLLc(k-1))
                else
                    y = (x - GLLc(k+1))/(GLLc(k)-GLLc(k+1))
                endif
            endif
        endif
    endif

    ! ################################################################
    return
end subroutine pol_force

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
