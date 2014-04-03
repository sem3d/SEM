!>
!!\file shape4.F90
!!\brief Contient la subroutine shape4.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module shape_lin

contains

    ! #########################################
    subroutine compute_normals(Tdomain,nf)

        use sdomain

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain
        integer, intent (IN) :: nf

        ! local variables
        real, dimension(0:1) :: normal_aux
        integer :: i, n_elem, ngll, nv0, nv1

        n_elem = Tdomain%sFace(nf)%Near_Element(0)
        ngll   = Tdomain%sFace(nf)%ngll
        nv0    = Tdomain%sFace(nf)%Near_Vertex(0)
        nv1    = Tdomain%sFace(nf)%Near_Vertex(1)
        allocate(Tdomain%sFace(nf)%Normal(0:1))
        allocate(Tdomain%sFace(nf)%Normal_Nodes(0:ngll-1,0:1))

        ! Computation of an unique Face normal
        call n_from_vertices(Tdomain, normal_aux, nv0, nv1)
        Tdomain%sFace(nf)%Normal = normal_aux
        ! Check if normal is outward...
        if (Tdomain%sFace(nf)%Which_Face(0) .GE. 2) then
            Tdomain%sFace(nf)%Normal(:) = - Tdomain%sFace(nf)%Normal(:)
        endif

        ! Computation of a field of normals for each Gauss points
        ! For shape4, elements are linears, and the normal for the
        ! Gauss points are the same than the Face Normal
        do i = 0,ngll-1
            Tdomain%sFace(nf)%Normal_Nodes(i,0) = Tdomain%sFace(nf)%Normal(0)
            Tdomain%sFace(nf)%Normal_Nodes(i,1) = Tdomain%sFace(nf)%Normal(1)
        end do

    end subroutine compute_normals


! #########################################
    subroutine compute_normals_elem(Tdomain, nelem)

        use sdomain

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain
        integer, intent (IN) :: nelem

        ! local variables
        integer              :: i, j, k, ngllx, ngllz, ngll, nglltot, nv0, nv1
        real, dimension(0:1) :: normal_aux

        ngllx  = Tdomain%specel(nelem)%ngllx
        ngllz  = Tdomain%specel(nelem)%ngllz
        nglltot= 2*(ngllx+ngllz)
        k = 0

        allocate(Tdomain%specel(nelem)%Normal_nodes(0:nglltot-1,0:1))

        do i=0,3
            nv0  = Tdomain%specel(nelem)%Near_Vertex(i)
            nv1  = Tdomain%specel(nelem)%Near_Vertex(mod(i+1,4))
            call n_from_vertices(Tdomain,normal_aux,nv0,nv1)
            if (i==0 .OR. i==2) then
                ngll = ngllx
            else
                ngll = ngllz
            endif
            do j=0,ngll-1
                Tdomain%specel(nelem)%Normal_Nodes(k,:) = normal_aux(:)
            enddo
            k = k+1
        enddo

    end subroutine compute_normals_elem


    ! #########################################
    !>
    !!\file shape4.F90
    !!\brief
    !!\version 1.0
    !!\date 01/04/2014
    !! This function computes a normal from 2 vertices
    !! Be carefull : this normal may not be outward !
    !<
    subroutine n_from_vertices(Tdomain,n_vec,nv0,nv1)

        use sDomain
        implicit none

        type (domain),target, intent (IN)    :: Tdomain
        real, dimension(0:1), intent (INOUT) :: n_vec
        integer, intent(IN)        :: nv0, nv1
        integer                    :: n0, n1
        real                       :: dx, dy, nx, ny, norm

        n0 = Tdomain%sVertex(nv0)%Glob_numbering
        n1 = Tdomain%sVertex(nv1)%Glob_numbering
        dx = Tdomain%coord_nodes(0,nv1) - Tdomain%coord_nodes(0,nv0)
        dy = Tdomain%coord_nodes(1,nv1) - Tdomain%coord_nodes(1,nv0)
        nx  = dy ; ny = -dx
        norm = sqrt(nx*nx + ny*ny)
        n_vec(0) = nx / norm
        n_vec(1) = ny / norm

    end subroutine n_from_vertices


    ! #########################################
    !>
    !!\file shape4.F90
    !!\brief
    !!\version 1.0
    !!\date 17/01/2014
    !! This subroutine computes the Jacobian of the transformation
    !! from a reference segment [-1,1] to the actual face nf. This
    !! tranformation is, in the linear case, just homothetic, therefore
    !! the Jacobian is the same for all the nodes of the face, and its
    !! value is equal to elongation coefficient.
    !<
    subroutine compute_Jacobian_1D(Tdomain,nf,Jac1D)

        use sdomain

        implicit none

        type(domain), intent (IN) :: Tdomain
        integer, intent (IN) :: nf
        real, intent (INOUT) :: Jac1D

        ! local variables
        integer :: nv0, nv1, n0, n1

        nv0= Tdomain%sFace(nf)%Near_Vertex(0)
        nv1= Tdomain%sFace(nf)%Near_Vertex(1)
        n0 = Tdomain%sVertex(nv0)%Glob_numbering
        n1 = Tdomain%sVertex(nv1)%Glob_numbering

        Jac1D = 0.5 * sqrt((Tdomain%Coord_Nodes(0,n0)-Tdomain%Coord_Nodes(0,n1))**2 &
            +(Tdomain%Coord_Nodes(1,n0)-Tdomain%Coord_Nodes(1,n1))**2 )

    end subroutine compute_Jacobian_1D

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
        integer :: nf
        real :: Jac1D

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

            if (Tdomain%specel(n)%Type_DG .NE. GALERKIN_CONT) then
                if (ngllx .NE. ngllz) STOP 'Case ngllx not equal to ngllz is not taken into account'
                allocate(Tdomain%specel(n)%Coeff_Integr_Faces(0:3,0:ngllx-1))
                do i=0,3
                    nf = Tdomain%specel(n)%Near_Face(i)
                    call compute_Jacobian_1D(Tdomain,nf,Jac1D)
                    Tdomain%specel(n)%Coeff_Integr_Faces(i,:) = Tdomain%sSubdomain(mat)%GLLwx(:) * Jac1D
                enddo
            endif
        enddo

        ! Compute Normals for Faces with DG :
        do nf=0,Tdomain%n_face-1
            if (Tdomain%sFace(nf)%Type_flux .NE. FLUX_NONE) then
                call compute_normals(Tdomain,nf)
            end if
        end do

        ! Compute Normals for Elements HDG :
        do n=0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%Type_DG .EQ. GALERKIN_HDG_RP) then
                call compute_normals_elem(Tdomain,n)
            end if
        end do

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
