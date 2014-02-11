!>
!!\file shape4.F90
!!\brief Contient la subroutine shape4.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Evalue les fonctions de forme.
!!  Ecriture des fichiers posPtg.*** et connecPtg.****
!! \param type(domain),target, intent (INOUT) Tdomain
!<

! #########################################
module shape_lin


contains

! #########################################

subroutine shape4(Tdomain)

    use sdomain
    use semdatafiles

    implicit none

    type(domain),target, intent (INOUT) :: Tdomain


    ! local variables

    integer :: i_aus,n, mat,ngllx,ngllz,i,j,ipoint
    real :: x0,x1,x2,x3,z0,z1,z2,z3,xi,eta,xp,zp, Jac
    real, dimension (0:1,0:1) :: LocInvGrad

    ! Modified by Gaetano Festa, 26/05/05
    !---------------------------------------------------------------------------------------------------------------
    ! Shape functions are derived from "Finite Elements and Approximations"
    ! by Zienkiewicz, O.C. and Morgan, K.
    ! John Wiley and Sons, 1983
    ! --------------------------------------------------------------------------------------------------------------

    allocate (Tdomain%GlobCoord(0:1,0:Tdomain%n_glob_points-1))

    do n = 0,Tdomain%n_elem - 1
        i_aus = Tdomain%specel(n)%Control_Nodes(0);  x0 = Tdomain%Coord_Nodes(0,i_aus);  z0 = Tdomain%Coord_Nodes(1,i_aus)
        i_aus = Tdomain%specel(n)%Control_Nodes(1);  x1 = Tdomain%Coord_Nodes(0,i_aus);  z1 = Tdomain%Coord_Nodes(1,i_aus)
        i_aus = Tdomain%specel(n)%Control_Nodes(2);  x2 = Tdomain%Coord_Nodes(0,i_aus);  z2 = Tdomain%Coord_Nodes(1,i_aus)
        i_aus = Tdomain%specel(n)%Control_Nodes(3);  x3 = Tdomain%Coord_Nodes(0,i_aus);  z3 = Tdomain%Coord_Nodes(1,i_aus)

        mat = Tdomain%specel(n)%mat_index
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz

        allocate (Tdomain%specel(n)%Jacob(0:ngllx-1,0:ngllz-1) )
        allocate (Tdomain%specel(n)%InvGrad(0:ngllx-1,0:ngllz-1,0:1,0:1) )

        do j = 0,ngllz - 1

            eta =   Tdomain%sSubdomain(mat)%GLLcz (j)

            do i = 0,ngllx - 1

                xi = Tdomain%sSubdomain(mat)%GLLcx (i)

                ipoint = Tdomain%specel(n)%Iglobnum(i,j)


                xp = 0.25 * ( x0 * (1.-xi)*(1.-eta) + x1 * (1.+xi)*(1.-eta) + x2 * (1.+xi)*(1.+eta) + x3 * (1.-xi)*(1.+eta) )
                zp = 0.25 * ( z0 * (1.-xi)*(1.-eta) + z1 * (1.+xi)*(1.-eta) + z2 * (1.+xi)*(1.+eta) + z3 * (1.-xi)*(1.+eta) )

                Tdomain%GlobCoord (0,ipoint) = xp;   Tdomain%GlobCoord (1,ipoint) = zp

                !         Computation of the derivative matrix, dx_(jj)/dxi_(ii)

                LocInvGrad(0,0) = 0.25 * ( (x1-x0) * (1-eta) + (x2-x3) * (1+eta) )
                LocInvGrad(1,0) = 0.25 * ( (x3-x0) * (1-xi) + (x2-x1) * (1+xi) )
                LocInvGrad(0,1) = 0.25 * ( (z1-z0) * (1-eta) + (z2-z3) * (1+eta) )
                LocInvGrad(1,1) = 0.25 * ( (z3-z0) * (1-xi) + (z2-z1) * (1+xi) )

                call invert2 (LocInvGrad, Jac)

                Tdomain%specel(n)%InvGrad (i,j,0:1,0:1)  = LocInvGrad (0:1,0:1)

                Tdomain%specel(n)%Jacob (i,j) = Jac
            enddo
        enddo

    enddo


    if (Tdomain%logicD%super_object_local_present) then
        do n = 0, Tdomain%n_fault-1
            call shape4_manage_super_object(Tdomain, n)
        enddo
    endif
    return
end subroutine shape4


subroutine shape4_manage_super_object(Tdomain, n)
    use sdomain
    implicit none
    type(domain),target, intent (INOUT) :: Tdomain
    integer, intent(IN) :: n
    integer :: i, j, i_aus
    real, dimension (:), pointer :: GLLc_face
    real, dimension (:,:), allocatable :: Store_normal
    real :: ds_local, normal_0, normal_1, normalization
    real :: x0, x1, z0, z1
    integer :: Face_up, mat, ngll, nv, nv2

            do  j = 0, Tdomain%sFault(n)%n_face-1
                ngll = Tdomain%sFault(n)%fFace(j)%ngll
                Face_up = Tdomain%sFault(n)%fFace(j)%Face_UP
                nv2 = Tdomain%sFace(Face_up)%Near_element(0)
                mat = Tdomain%specel(nv2)%mat_index
                if (Tdomain%specel(nv2)%near_face(0) == Face_up .or. Tdomain%specel(nv2)%near_face(2) == face_up) then
                    GLLc_face => Tdomain%sSubdomain(mat)%GLLcx
                else
                    GLLc_face => Tdomain%sSubdomain(mat)%GLLcz
                endif

                if (Tdomain%specel(nv2)%near_face(0) == Face_up) then
                    i_aus = Tdomain%specel(nv2)%Control_Nodes(0)
                    x0 = Tdomain%Coord_Nodes(0,i_aus)
                    z0 = Tdomain%Coord_Nodes(1,i_aus)
                    i_aus = Tdomain%specel(nv2)%Control_Nodes(1)
                    x1 = Tdomain%Coord_Nodes(0,i_aus)
                    z1 = Tdomain%Coord_Nodes(1,i_aus)
                else if (Tdomain%specel(nv2)%near_face(1) == Face_up) then
                    i_aus = Tdomain%specel(nv2)%Control_Nodes(1)
                    x0 = Tdomain%Coord_Nodes(0,i_aus)
                    z0 = Tdomain%Coord_Nodes(1,i_aus)
                    i_aus = Tdomain%specel(nv2)%Control_Nodes(2)
                    x1 = Tdomain%Coord_Nodes(0,i_aus)
                    z1 = Tdomain%Coord_Nodes(1,i_aus)
                else if (Tdomain%specel(nv2)%near_face(2) == Face_up) then
                    i_aus = Tdomain%specel(nv2)%Control_Nodes(2)
                    x0 = Tdomain%Coord_Nodes(0,i_aus)
                    z0 = Tdomain%Coord_Nodes(1,i_aus)
                    i_aus = Tdomain%specel(nv2)%Control_Nodes(3)
                    x1 = Tdomain%Coord_Nodes(0,i_aus)
                    z1 = Tdomain%Coord_Nodes(1,i_aus)
                else if (Tdomain%specel(nv2)%near_face(3) == Face_up) then
                    i_aus = Tdomain%specel(nv2)%Control_Nodes(3)
                    x0 = Tdomain%Coord_Nodes(0,i_aus)
                    z0 = Tdomain%Coord_Nodes(1,i_aus)
                    i_aus = Tdomain%specel(nv2)%Control_Nodes(0)
                    x1 = Tdomain%Coord_Nodes(0,i_aus)
                    z1 = Tdomain%Coord_Nodes(1,i_aus)
                else
                    stop "Inconsistency with Face_up"
                endif

                allocate (Tdomain%sFault(n)%fFace(j)%ds(0:ngll-1))
                allocate (Tdomain%sFault(n)%fFace(j)%normal(0:ngll-1,0:1))
                allocate (Tdomain%sFault(n)%fFace(j)%distance(0:ngll-1))

                ds_local =  (x1-x0)**2 + (z1-z0)**2
                ds_local =  Sqrt(ds_local)
                normal_0 = (z1 - z0)/ds_local; normal_1 = (x0-x1)/ds_local

                do  i = 0, ngll-1
                    Tdomain%sFault(n)%fFace(j)%distance(i) = ds_local * (1+GLLc_face(i))
                enddo

                Tdomain%sFault(n)%fFace(j)%X_Vertex(0) = x0; Tdomain%sFault(n)%fFace(j)%X_Vertex(1) = x1
                Tdomain%sFault(n)%fFace(j)%Z_Vertex(0) = z0; Tdomain%sFault(n)%fFace(j)%Z_Vertex(1) = z1
                Tdomain%sFault(n)%fFace(j)%ds(:) = 0.5 * ds_local
                Tdomain%sFault(n)%fFace(j)%normal(0:ngll-1,0) = normal_0; Tdomain%sFault(n)%fFace(j)%normal(0:ngll-1,1) = normal_1
            enddo

            do j = 0, Tdomain%sFault(n)%n_vertex-1
                Tdomain%sFault(n)%fVertex(j)%normal(0:1) = 0
            enddo

            do j = 0, Tdomain%sFault(n)%n_face-1
                ngll = Tdomain%sFault(n)%fFace(j)%ngll
                nv = Tdomain%sFault(n)%fFace(j)%Face_To_Vertex(0)
                Tdomain%sFault(n)%fVertex(nv)%Normal(0:1) = Tdomain%sFault(n)%fVertex(nv)%Normal(0:1) + &
                    Tdomain%sFault(n)%fFace(j)%normal(0,0:1)
                nv = Tdomain%sFault(n)%fFace(j)%Face_To_Vertex(1)
                Tdomain%sFault(n)%fVertex(nv)%Normal(0:1) = Tdomain%sFault(n)%fVertex(nv)%Normal(0:1) + &
                    Tdomain%sFault(n)%fFace(j)%normal(ngll-1,0:1)
            enddo

            do j = 0, Tdomain%sFault(n)%n_vertex-1
                normalization = Tdomain%sFault(n)%fVertex(j)%normal(0)**2+Tdomain%sFault(n)%fVertex(j)%normal(1)**2
                normalization = sqrt(normalization)
                Tdomain%sFault(n)%fVertex(j)%normal(:) = Tdomain%sFault(n)%fVertex(j)%normal(:)/normalization
            enddo;

            do j = 0, Tdomain%sFault(n)%n_face-1
                ngll = Tdomain%sFault(n)%fFace(j)%ngll
                allocate (Store_normal(1:ngll-2,0:1))
                Store_normal (1:ngll-2,0:1) = Tdomain%sFault(n)%fFace(j)%normal(1:ngll-2,0:1)
                deallocate(Tdomain%sFault(n)%fFace(j)%normal)
                allocate(Tdomain%sFault(n)%fFace(j)%normal(1:ngll-2,0:1))
                Tdomain%sFault(n)%fFace(j)%normal(1:ngll-2,0:1) = Store_normal (1:ngll-2,0:1)
                deallocate (Store_normal)
            enddo

            ! Valid only if points are ordered
            ds_local = 0
            do j =  0, Tdomain%sFault(n)%n_face-1
                ngll = Tdomain%sFault(n)%fFace(j)%ngll
                do i= 0, ngll-1
                    Tdomain%sFault(n)%fFace(j)%distance(i) = ds_local +  Tdomain%sFault(n)%fFace(j)%distance(i)
                enddo
                ds_local = ds_local + (Tdomain%sFault(n)%fFace(j)%distance(ngll-1) - Tdomain%sFault(n)%fFace(j)%distance(0))
                nv = Tdomain%sFault(n)%fFace(j)%Face_to_vertex(0)
                Tdomain%sFault(n)%fvertex(nv)%distance = Tdomain%sFault(n)%fFace(j)%distance(0)
                nv = Tdomain%sFault(n)%fFace(j)%Face_to_vertex(1)
                Tdomain%sFault(n)%fvertex(nv)%distance = Tdomain%sFault(n)%fFace(j)%distance(ngll-1)
                allocate (Store_normal(0:ngll-1,0:0))
                Store_normal(:,0) =  Tdomain%sFault(n)%fFace(j)%distance
                deallocate ( Tdomain%sFault(n)%fFace(j)%distance)
                allocate ( Tdomain%sFault(n)%fFace(j)%distance(1:ngll-2))
                Tdomain%sFault(n)%fFace(j)%distance(1:ngll-2)= Store_normal (1:ngll-2,0)
                deallocate (Store_normal)
            enddo

end subroutine shape4_manage_super_object

subroutine shape4_local_coord(xc, zc, x, z, xi1, eta1, inosol)
    implicit none
    real, dimension (0:3), intent(in) :: xc,zc
    real, intent(in) :: x, z
    real, intent(out) :: xi1, eta1
    logical, intent(out) :: inosol
    !
    integer :: i
    real :: a1, b1, c1, d1
    real :: a2, b2, c2, d2
    real :: alpha, beta, gamm, delta

    a1 =  4 * x - xc(0) - xc(1) - xc(2) - xc(3)
    b1 =  xc(0) - xc(1) + xc(3) - xc(2)
    c1 =  xc(0) + xc(1) - xc(2) - xc(3)
    d1 = -xc(0) + xc(1) + xc(3) - xc(2)
    a2 =  4 * z - zc(0) - zc(1) - zc(2) - zc(3)
    b2 =  zc(0) - zc(1) + zc(3) - zc(2)
    c2 =  zc(0) + zc(1) - zc(2) - zc(3)
    d2 = -zc(0) + zc(1) + zc(3) - zc(2)
    alpha = c1*d2 - d1*c2  ; beta = a1*d2 - b1*c2 + c1*b2 - d1*a2; gamm = a1*b2 - a2*b1
    if (abs(alpha)<1e-7 ) then
        eta1 = -gamm/beta
        if (d2 == 0 .and. b2==0) then
            xi1 = -(a1 + c1*eta1)/(b1+d1*eta1)
        else
            xi1 = -(a2 + c2*eta1)/(b2+d2*eta1)
        endif
        inosol = xi1 <=1 .and. xi1>=-1 .and. eta1>=-1 .and. eta1<=1
        inosol =.not. inosol
    else
        delta = beta**2 - 4* alpha*gamm
        if (delta < 0) then
            write (*,*)  "No solution for the location"
            write (*,*) " Return to continue, and Ctrl C to quit"
            stop
        endif
        eta1 = 0.5 * (- beta + sqrt (delta) )/ alpha
        inosol = .true.
        if (eta1 <= 1 .and. eta1 >=-1) then
            xi1 = -(a2 + c2*eta1)/(b2+d2*eta1)
            if (xi1 <=1 .and. xi1 >= -1) inosol = .false.
        endif
        if (inosol) then
            eta1 =  0.5 * (- beta - sqrt (delta) )/ alpha
            if (eta1 <= 1 .and. eta1 >=-1) then
                xi1 = -(a2 + c2*eta1)/(b2+d2*eta1)
                if (xi1 <=1 .and. xi1 >= -1) inosol = .false.
            endif
        endif
    endif
    if (inosol) then
        write (*,*)  "No solution found for coordinates    ",x, z
        write (*,*)  "Within element :"
        do i=0,3
            write(*,*) xc(i), zc(i)
        end do
        stop
    endif

end subroutine shape4_local_coord

end module shape_lin
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
