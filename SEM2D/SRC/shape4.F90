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

! #########################################

subroutine compute_normals(Tdomain,nf)

    use sdomain

    implicit none

    type(domain),target, intent (INOUT) :: Tdomain
    integer, intent (IN) :: nf

    ! local variables
    integer :: i, n_elem, ngll, nv0, nv1
    real    :: n1, n2, nx, ny, norm

    n_elem = Tdomain%sFace(nf)%Near_Element(0)
    ngll   = Tdomain%sFace(nf)%ngll
    nv0    = Tdomain%sFace(nf)%Near_Vertex(0)
    nv1    = Tdomain%sFace(nf)%Near_Vertex(1)
    allocate(Tdomain%sFace(nf)%Normal(0:1))
    allocate(Tdomain%sFace(nf)%Normal_Nodes(0:ngll-1,0:1))

    ! Computation of an unique Face normal
    nv0 = Tdomain%sVertex(nv0)%Glob_numbering
    nv1 = Tdomain%sVertex(nv1)%Glob_numbering
    n1  = Tdomain%coord_nodes(0,nv1) - Tdomain%coord_nodes(0,nv0)
    n2  = Tdomain%coord_nodes(1,nv1) - Tdomain%coord_nodes(1,nv0)
    nx  = n2 ; ny = -n1
    if (Tdomain%sFace(nf)%Which_Face(0) .GE. 2) then
        nx = -nx ; ny = -ny
    endif
    norm = sqrt(nx*nx + ny*ny)
    Tdomain%sFace(nf)%Normal(0) = nx / norm
    Tdomain%sFace(nf)%Normal(1) = ny / norm

    ! Computation of a field of normals for each Gauss points
    ! For shape4, elements are linears, and the normal for the
    ! Gauss points are the same than the Face Normal
    do i = 0,ngll-1
       Tdomain%sFace(nf)%Normal_Nodes(i,0) = Tdomain%sFace(nf)%Normal(0)
       Tdomain%sFace(nf)%Normal_Nodes(i,1) = Tdomain%sFace(nf)%Normal(1)
    end do

end subroutine compute_normals


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

!    integer :: i_aus,n, mat,ngllx,ngllz,i,j,ipoint, ngll, nv, Face_UP, nv2, nf
!    real :: x0,x1,x2,x3,z0,z1,z2,z3,xi,eta,xp,zp, Jac, ds_local, Jac1D
!    real :: normal_0, normal_1, normalization
    integer :: i_aus,n, mat,ngllx,ngllz,i,j,ipoint, nf
    real :: x0,x1,x2,x3,z0,z1,z2,z3,xi,eta,xp,zp, Jac
    real :: Jac1D
    real :: face_len

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

        if (Tdomain%specel(n)%Type_DG .NE. 2) then
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
       if (Tdomain%sFace(nf)%Type_flux .NE. 0) then
          call compute_normals(Tdomain,nf)
       end if
    end do

    if (Tdomain%logicD%super_object_local_present) then
        do n = 0, Tdomain%n_fault-1
            call manage_super_object(Tdomain, n)
        enddo
    endif
    return
end subroutine shape4


subroutine manage_super_object(Tdomain, n)
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

end subroutine manage_super_object

end module shape_lin
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
