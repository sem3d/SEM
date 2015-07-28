!>
!!\file shape8.F90
!!\brief Contient la subroutine shape8.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Evalue les fonctions de forme.
!!
!! \param type(domain),target, intent (INOUT) Tdomain
!<

module shape_quad
    use constants
    implicit none
contains

    subroutine shape8(Tdomain)

        use sdomain
        use semdatafiles
        use shape_lin ! compute_normals as quad4

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain


        ! local variables

        integer :: i_aus,n, mat,ngllx,ngllz,i,j,ipoint
        real, dimension(0:7) :: xc, zc
        real :: xi,eta,xp,zp, Jac
        real, dimension (0:1,0:1) :: LocInvGrad


        ! Modified by Gaetano Festa, 27/05/05
        !---------------------------------------------------------------------------
        ! Shape functions are derived from "Finite Elements and Approximations"
        ! by Zienkiewicz, O.C. and Morgan, K.
        ! John Wiley and Sons, 1983
        ! --------------------------------------------------------------------------

        allocate (Tdomain%GlobCoord(0:1,0:Tdomain%n_glob_points-1))

        do n = 0,Tdomain%n_elem - 1
            do i=0,7
                i_aus = Tdomain%specel(n)%Control_Nodes(i)
                xc(i) = Tdomain%Coord_Nodes(0,i_aus)
                zc(i) = Tdomain%Coord_Nodes(1,i_aus)
            end do

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

                    call shape8_global_coord(xc, zc, xi, eta, xp, zp)

                    Tdomain%GlobCoord (0,ipoint) = xp;   Tdomain%GlobCoord (1,ipoint) = zp

                    !     Computation of the derivative matrix, dx_(jj)/dxi_(ii) : more precisely we have :
                    !  LocInvGrad(0,0) = dx/dxi  ;  LocInvGrad(1,0) = dx/deta
                    !  LocInvGrad(0,1) = dz/dxi  ;  LocInvGrad(1,1) = dz/deta

                    LocInvGrad(0,0) = 0.25 * (xc(0) *(1-eta)*(2*xi+eta) &
                        +xc(1) *(1-eta)*(2*xi-eta) &
                        +xc(2) *(1+eta)*(2*xi+eta) &
                        +xc(3) *(1+eta)*(2*xi-eta)) &
                        -xc(4) * xi*(1-eta) &
                        -xc(6) * xi*(1+eta) &
                        +0.5* (xc(5)-xc(7))* (1-eta**2)
                    LocInvGrad(1,0) = 0.25 * (xc(0) *(1-xi)*(2*eta+xi) &
                        -xc(1) *(1+xi)*(xi-2*eta) &
                        +xc(2) *(1+xi)*(2*eta+xi) &
                        -xc(3) *(1-xi)*(xi-2*eta)) &
                        -xc(5) *eta*(1+xi) &
                        -xc(7) *eta*(1-xi) &
                        +0.5* (xc(6)-xc(4))* (1-xi**2)
                    LocInvGrad(0,1) = 0.25 * (zc(0) *(1-eta)*(2*xi+eta) &
                        +zc(1) *(1-eta)*(2*xi-eta) &
                        +zc(2) *(1+eta)*(2*xi+eta) &
                        +zc(3) *(1+eta)*(2*xi-eta)) &
                        -zc(4) * xi*(1-eta) &
                        -zc(6) *xi *(1+eta) &
                        +0.5* (zc(5)-zc(7))* (1-eta**2)
                    LocInvGrad(1,1) = 0.25 * (zc(0) *(1-xi)*(2*eta+xi) &
                        -zc(1) *(1+xi)*(xi-2*eta) &
                        +zc(2) *(1+xi)*(2*eta+xi) &
                        -zc(3) *(1-xi)*(xi-2*eta)) &
                        -zc(5) *eta*(1+xi) &
                        -zc(7) *eta *(1-xi) &
                        +0.5* (zc(6)-zc(4))* (1-xi**2)

                    !  Creation of the normal if the node is on an element side

                    if ((i.EQ.0) .OR. (j.EQ.0) .OR. (i.EQ.ngllx-1) .OR. (j.EQ.ngllz-1)) then
                        ! call buildNormal(Tdomain,LocInvGrad,n,i,j)
                    end if

                    call invert2 (LocInvGrad, Jac )

                    !  Computation of the local Jacobian and inversion of the Jacobian Matrix
                    Tdomain%specel(n)%InvGrad (i,j,0:1,0:1)  = LocInvGrad (0:1,0:1)

                    Tdomain%specel(n)%Jacob (i,j) = Jac
                enddo
            enddo
        enddo

        if (Tdomain%logicD%super_object_local_present) then
            do n = 0, Tdomain%n_fault-1
                call shape8_manage_super_object(Tdomain, n)
            end do
        endif

        ! Testing if the normals are defined for all the faces
        do n=0,Tdomain%n_face-1
            call compute_normals(Tdomain,n)
            ngllx =  Tdomain%sFace(n)%ngll
            do i=0,ngllx-1
                if((Tdomain%sFace(n)%Normal_Nodes(i,0).EQ.0.) .AND. (Tdomain%sFace(n)%Normal_Nodes(i,1).EQ.0.)) then
                    STOP "The Normal of the faces are not Computed Properly"
                endif
            end do
        end do

        return
    end subroutine shape8


    subroutine shape8_manage_super_object(Tdomain, n)
        use sdomain
        implicit none
        type(domain),target, intent (INOUT) :: Tdomain
        integer, intent(IN) :: n
        integer :: i, j, i_aus
        real, dimension (:), pointer :: GLLc_face
        real, dimension (:,:), allocatable :: Store_normal
        real :: ds_local, normal_0, normal_1, normalization
        real :: x0, x1, x2, z0, z1, z2, xp, zp
        integer :: Face_up, mat, ngll, nv, nv2
        real :: zin,xin,loc_distance


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
                i_aus = Tdomain%specel(nv2)%Control_Nodes(0);
                x0 = Tdomain%Coord_Nodes(0,i_aus);  z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(1);
                x1 = Tdomain%Coord_Nodes(0,i_aus);  z1 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(4);
                x2 = Tdomain%Coord_Nodes(0,i_aus);  z2 = Tdomain%Coord_Nodes(1,i_aus)
            else if (Tdomain%specel(nv2)%near_face(1) == Face_up) then
                i_aus = Tdomain%specel(nv2)%Control_Nodes(1);
                x0 = Tdomain%Coord_Nodes(0,i_aus);  z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(2);
                x1 = Tdomain%Coord_Nodes(0,i_aus);  z1 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(5);
                x2 = Tdomain%Coord_Nodes(0,i_aus);  z2 = Tdomain%Coord_Nodes(1,i_aus)
            else if (Tdomain%specel(nv2)%near_face(2) == Face_up) then
                i_aus = Tdomain%specel(nv2)%Control_Nodes(2);
                x0 = Tdomain%Coord_Nodes(0,i_aus);  z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(3);
                x1 = Tdomain%Coord_Nodes(0,i_aus);  z1 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(6);
                x2 = Tdomain%Coord_Nodes(0,i_aus);  z2 = Tdomain%Coord_Nodes(1,i_aus)
            else
                i_aus = Tdomain%specel(nv2)%Control_Nodes(3);
                x0 = Tdomain%Coord_Nodes(0,i_aus);  z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(0);
                x1 = Tdomain%Coord_Nodes(0,i_aus);  z1 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(7);
                x2 = Tdomain%Coord_Nodes(0,i_aus);  z2 = Tdomain%Coord_Nodes(1,i_aus)
            endif

            allocate (Tdomain%sFault(n)%fFace(j)%ds(0:ngll-1))
            allocate (Tdomain%sFault(n)%fFace(j)%normal(0:ngll-1,0:1))
            allocate (Tdomain%sFault(n)%fFace(j)%distance(0:ngll-1))

            loc_distance = 0.
            xin = x0; zin = z0

            do i = 0, ngll-1
                xp = x1 * (GLLc_face(i)+0.5) + x0 * (GLLc_face(i)-0.5) - 2*x2* GLLc_face(i)
                zp  = z1 * (GLLc_face(i)+0.5) + z0 * (GLLc_face(i)-0.5) - 2*z2* GLLc_face(i)

                ds_local = xp**2 + zp**2
                ds_local = sqrt(ds_local)
                normal_0 = zp/ds_local; normal_1 = -xp/ds_local
                Tdomain%sFault(n)%fFace(j)%ds(i) = ds_local
                Tdomain%sFault(n)%fFace(j)%normal(i,0) = normal_0; Tdomain%sFault(n)%fFace(j)%normal(i,1) = normal_1

                Tdomain%sFault(n)%fFace(j)%distance(i) = loc_distance
                if (i < ngll-1) then
                    xp = 0.5 * GLLc_face(i+1)*(GLLc_face(i+1)+1) *x1 +  0.5 * GLLc_face(i+1)*(-GLLc_face(i+1)+1) *x0 + &
                        x2 * (1 -  GLLc_face(i+1)**2)
                    zp = 0.5 * GLLc_face(i+1)*(GLLc_face(i+1)+1) *z1 +  0.5 * GLLc_face(i+1)*(-GLLc_face(i+1)+1) *z0 + &
                        z2 * (1 -  GLLc_face(i+1)**2)
                    loc_distance = loc_distance + sqrt((xp-xin)**2 + (zp-zin)**2)
                    xin = xp; zin = zp
                endif
            enddo
            Tdomain%sFault(n)%fFace(j)%X_Vertex(0) = x0; Tdomain%sFault(n)%fFace(j)%X_Vertex(1) = x1
            Tdomain%sFault(n)%fFace(j)%Z_Vertex(0) = z0; Tdomain%sFault(n)%fFace(j)%Z_Vertex(1) = z1
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
        enddo

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
            ds_local = Tdomain%sFault(n)%fFace(j)%distance(ngll-1)
            nv = Tdomain%sFault(n)%fFace(j)%Face_to_vertex(0)
            Tdomain%sFault(n)%fvertex(nv)%distance = Tdomain%sFault(n)%fFace(j)%distance(0)
            nv = Tdomain%sFault(n)%fFace(j)%Face_to_vertex(1)
            Tdomain%sFault(n)%fvertex(nv)%distance = Tdomain%sFault(n)%fFace(j)%distance(ngll-1)
            allocate (Store_normal(0:ngll-1,0:0))
            Store_normal(:,0) =  Tdomain%sFault(n)%fFace(j)%distance
            allocate ( Tdomain%sFault(n)%fFace(j)%distance(1:ngll-2))
            Tdomain%sFault(n)%fFace(j)%distance(1:ngll-2)= Store_normal (1:ngll-2,0)
            deallocate (Store_normal)
        enddo

        return
    end subroutine shape8_manage_super_object

!##########################################

    subroutine buildNormal(Tdomain,LocInvGrad,n_elem,i,j)

        use sdomain

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain
        real, dimension (0:1,0:1), intent (IN) :: LocInvGrad
        integer, intent (IN) :: n_elem, i, j

        ! local variables
        integer :: ngllx, ngllz, nf
        real    :: tx, tz, nx, nz, norm_n

        ngllx = Tdomain%specel(n_elem)%ngllx
        ngllz = Tdomain%specel(n_elem)%ngllz

        if(i.EQ.0) then ! local face 3
            nf = Tdomain%specel(n_elem)%Near_Face(3)
            if(Tdomain%sface(nf)%Near_Element(0) .EQ. n_elem) then
                ! Tangent from derivatives along eta
                tx = LocInvGrad(1,0) ; tz = LocInvGrad(1,1)
                ! Normal build from tangent
                nx = -tz ; nz = tx
                norm_n = sqrt(nx**2 + nz**2)
                Tdomain%sFace(nf)%Normal_Nodes(j,0) = nx / norm_n
                Tdomain%sFace(nf)%Normal_Nodes(j,1) = nz / norm_n
            end if
        else if (i.EQ.ngllx-1) then ! local face 1
            nf = Tdomain%specel(n_elem)%Near_Face(1)
            if(Tdomain%sface(nf)%Near_Element(0) .EQ. n_elem) then
                ! Tangent from derivatives along eta
                tx = LocInvGrad(1,0) ; tz = LocInvGrad(1,1)
                ! Normal build from tangent
                nx = tz ; nz = -tx
                norm_n = sqrt(nx**2 + nz**2)
                Tdomain%sFace(nf)%Normal_Nodes(j,0) = nx / norm_n
                Tdomain%sFace(nf)%Normal_Nodes(j,1) = nz / norm_n
            end if
        end if

        if(j.EQ.0) then ! local face 0
            nf = Tdomain%specel(n_elem)%Near_Face(0)
            if(Tdomain%sface(nf)%Near_Element(0) .EQ. n_elem) then
                ! Tangent from derivatives along eta
                tx = LocInvGrad(0,0) ; tz = LocInvGrad(0,1)
                ! Normal build from tangent
                nx = tz ; nz = -tx
                norm_n = sqrt(nx**2 + nz**2)
                Tdomain%sFace(nf)%Normal_Nodes(i,0) = nx / norm_n
                Tdomain%sFace(nf)%Normal_Nodes(i,1) = nz / norm_n
            end if
        else if (j.EQ.ngllz-1) then ! local face 2
            nf = Tdomain%specel(n_elem)%Near_Face(2)
            if(Tdomain%sface(nf)%Near_Element(0) .EQ. n_elem) then
                ! Tangent from derivatives along eta
                tx = LocInvGrad(0,0) ; tz = LocInvGrad(0,1)
                ! Normal build from tangent
                nx = -tz ; nz = tx
                norm_n = sqrt(nx**2 + nz**2)
                Tdomain%sFace(nf)%Normal_Nodes(i,0) = nx / norm_n
                Tdomain%sFace(nf)%Normal_Nodes(i,1) = nz / norm_n
            end if
        end if

        return
    end subroutine buildNormal

    subroutine shape8_global_coord(xq, zq, xi, eta, x, z)
        real, dimension(0:7), intent(in) :: xq, zq
        real, intent(in) :: xi, eta
        real, intent(out) :: x, z

        x = 0.25 * ( -xq(0) * (1.-xi)*(1.-eta)*(1+xi+eta) &
            -xq(1) * (1.+xi)*(1.-eta)*(1-xi+eta) &
            -xq(2) * (1.+xi)*(1.+eta)*(1-xi-eta) &
            -xq(3) * (1.-xi)*(1.+eta)*(1+xi-eta) ) &
            + 0.5 * (  xq(4) * (1.-xi**2)*(1.-eta) + &
            xq(5) * (1.+xi)*(1.-eta**2) + &
            xq(6) * (1.-xi**2)*(1.+eta) + &
            xq(7) * (1.-xi)*(1.-eta**2) )
        z = 0.25 * ( -zq(0) * (1.-xi)*(1.-eta)*(1+xi+eta) &
            -zq(1) * (1.+xi)*(1.-eta)*(1-xi+eta) &
            -zq(2) * (1.+xi)*(1.+eta)*(1-xi-eta) &
            -zq(3) * (1.-xi)*(1.+eta)*(1+xi-eta) ) &
            + 0.5 * (  zq(4) * (1.-xi**2)*(1.-eta) + &
            zq(5) * (1.+xi)*(1.-eta**2) + &
            zq(6) * (1.-xi**2)*(1.+eta) + &
            zq(7) * (1.-xi)*(1.-eta**2) )

    end subroutine shape8_global_coord

    subroutine shape8_local_coord(xq, zq, x, z, xi1, eta1, inosol)
        implicit none
        real, dimension (0:7), intent(in) :: xq, zq
        real, intent(in) :: x, z
        real, intent(out) :: xi1, eta1
        logical, intent(out) :: inosol
        !
        real, dimension(0:7) :: xc, zc
        integer, parameter :: nimax = 50, njmax = 50
        real, parameter :: dximax = 2./nimax,  detamax = 2./njmax
        real :: xi2, eta2
        logical :: inner
        integer :: i, j

        inner = .false.
        do12_jmax : do j = 0,njmax-1
            do i = 0, nimax-1
                xi1 = i*dximax -1; xi2 = xi1 + dximax
                eta1 = j * detamax-1 ; eta2 = eta1 + detamax
                xi1 = xi1 -dximax/nimax; xi2 = xi2 + dximax/nimax;
                eta1 = eta1 - detamax/njmax; eta2 = eta2 + detamax/njmax

                call shape8_global_coord(xq, zq, xi1, eta1, xc(0), zc(0))
                call shape8_global_coord(xq, zq, xi2, eta1, xc(1), zc(1))
                call shape8_global_coord(xq, zq, xi2, eta2, xc(2), zc(2))
                call shape8_global_coord(xq, zq, xi1, eta2, xc(3), zc(3))

                call verify_in_quad (Xc, Zc, x, z, inner)
                if (inner) exit do12_jmax
            enddo
        enddo do12_jmax

        inosol = .not. inner
        xi1 = 0.5*(xi1+xi2)
        eta1 = 0.5*(eta1+eta2)
    end subroutine shape8_local_coord

end module shape_quad
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
