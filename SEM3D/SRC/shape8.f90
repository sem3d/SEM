!>
!!\file shape8.f90
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

#ifdef MKA3D
#ifdef USE_HDF5
subroutine save_post_data(Tdomain,rg)
    use sdomain
    use semdatafiles
    use HDF5
    use sem_hdf5, only : create_dset, init_hdf5, create_dset_2d
    implicit none
    type(domain),target, intent (INOUT) :: Tdomain
    integer, intent(IN) :: rg
    character (len=MAX_FILE_SIZE) :: fnamef
    integer :: hdferr, nglobnum, ig
    integer :: n, i, j, k, ngllx, nglly, ngllz
    integer(HID_T) :: fid, dset_id
    integer, dimension(3,0:Tdomain%n_elem-1) :: ngll
    integer, dimension(:), allocatable :: iglobnum
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2

    write(*,*) "Creation donnees posttraitement"
    call semname_post_data(rg,fnamef)
    call init_hdf5()
    call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

    ! Sauvegarde du nombre de points de gauss par maille
    call create_dset_2d(fid, "NGLL", H5T_STD_I16LE, 3, Tdomain%n_elem, dset_id)
    nglobnum = 0
    do n = 0,Tdomain%n_elem - 1
        ngll(1,n) = Tdomain%specel(n)%ngllx
        ngll(2,n) = Tdomain%specel(n)%nglly
        ngll(3,n) = Tdomain%specel(n)%ngllz
        nglobnum = nglobnum + ngll(1,n)*ngll(2,n)*ngll(3,n)
    end do
    dim2(1) = 3
    dim2(2) = Tdomain%n_elem
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ngll, dim2, hdferr)
    call h5dclose_f(dset_id, hdferr)

    ! Sauvegarde de la numerotation globale des points de gauss par element
    call create_dset(fid, "IGLOBNUM", H5T_STD_I32LE, nglobnum, dset_id)
    allocate (iglobnum(nglobnum))
    ig = 1
    do n = 0,Tdomain%n_elem - 1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        do k = 0,ngllz - 1
            do j = 0,nglly - 1
                do i = 0,ngllx - 1
                    iglobnum(ig) = Tdomain%specel(n)%Iglobnum(i,j,k)
                    ig = ig + 1
                end do
            end do
        end do
    end do
    dim1(1) = nglobnum
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, iglobnum, dim1, hdferr)
    deallocate(iglobnum)
    call h5dclose_f(dset_id, hdferr)

    ! Sauvegarde des positions des points de gauss
    call create_dset_2d(fid, "NODES", H5T_IEEE_F64LE, 3, Tdomain%n_glob_points, dset_id)
    dim2(1) = 3
    dim2(2) = Tdomain%n_glob_points
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Tdomain%GlobCoord, dim2, hdferr)
    call h5dclose_f(dset_id, hdferr)
    call h5fclose_f(fid, hdferr)
    write(*,*) "fin Creation donnees posttraitement"
end subroutine save_post_data
#else
subroutine save_post_data(Tdomain,rg)
    use sdomain
    implicit none
    type(domain),target, intent (INOUT) :: Tdomain
    integer, intent(IN) :: rg, n, i, j, k, ipoint
    character (len=MAX_FILE_SIZE) :: fnamef

    call semname_connecptg(rg,fnamef)
    open (89,file=fnamef,form="formatted",status="unknown")
    do n = 0,Tdomain%n_elem - 1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        write(89,'(I8,3(X,I4))') n, ngllx, nglly, ngllz

        do k = 0,ngllz - 1
            do j = 0,nglly - 1
                do i = 0,ngllx - 1
                    ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)
                    write(89,*) ipoint
                end do
            end do
        end do
    end do
    close(89)

    call semname_posptg(rg,fnamef)
    open (88,file=fnamef,form="formatted",status="unknown")
    do i=0,Tdomain%n_glob_points-1
        write(88,'(I9,1X,3(1pe15.8,1X))') i, (Tdomain%GlobCoord (j,i),j=0,2)
    enddo
    close(88)
end subroutine save_post_data
#endif
#else
subroutine save_post_data(Tdomain,rg)
    use sdomain
    implicit none
    type(domain),target, intent (INOUT) :: Tdomain
    integer, intent(IN) :: rg
end subroutine save_post_data
#endif

subroutine shape8(Tdomain,rg)

    ! Modified by Paul Cupillard 08/12/2005
    use sdomain
    use semdatafiles

    implicit none

    type(domain),target, intent (INOUT) :: Tdomain
    integer, intent(IN) :: rg
    integer :: n, i_aus, ngll,ngllx,nglly,ngllz,ngll1,ngll2, mat, i,j,k ,ipoint, nf,ne,nv,n_aus,n_up,n_down
    real :: x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x7, y7, z7
    real :: xi,eta,zeta, xp,yp,zp, Jac,norm
    real, dimension (0:2,0:1) :: LocInvGradS
    real, dimension (0:2,0:2) :: LocInvGrad
    character (len=MAX_FILE_SIZE) :: fnamef

    !---------------------------------------------------------------------------------------------------------------
    ! Shape functions are derived from "The finite element displayed"
    ! by Dhatt, G. and Touzot, G.
    ! John Wiley and Sons, 1984
    ! --------------------------------------------------------------------------------------------------------------
    allocate (Tdomain%GlobCoord(0:2,0:Tdomain%n_glob_points-1))

    do n = 0,Tdomain%n_elem - 1

        i_aus = Tdomain%specel(n)%Control_Nodes(0)
        x0 = Tdomain%Coord_Nodes(0,i_aus)
        y0 = Tdomain%Coord_Nodes(1,i_aus)
        z0 = Tdomain%Coord_Nodes(2,i_aus)

        i_aus = Tdomain%specel(n)%Control_Nodes(1)
        x1 = Tdomain%Coord_Nodes(0,i_aus)
        y1 = Tdomain%Coord_Nodes(1,i_aus)
        z1 = Tdomain%Coord_Nodes(2,i_aus)

        i_aus = Tdomain%specel(n)%Control_Nodes(2)
        x2 = Tdomain%Coord_Nodes(0,i_aus)
        y2 = Tdomain%Coord_Nodes(1,i_aus)
        z2 = Tdomain%Coord_Nodes(2,i_aus)

        i_aus = Tdomain%specel(n)%Control_Nodes(3)
        x3 = Tdomain%Coord_Nodes(0,i_aus)
        y3 = Tdomain%Coord_Nodes(1,i_aus)
        z3 = Tdomain%Coord_Nodes(2,i_aus)

        i_aus = Tdomain%specel(n)%Control_Nodes(4)
        x4 = Tdomain%Coord_Nodes(0,i_aus)
        y4 = Tdomain%Coord_Nodes(1,i_aus)
        z4 = Tdomain%Coord_Nodes(2,i_aus)

        i_aus = Tdomain%specel(n)%Control_Nodes(5)
        x5 = Tdomain%Coord_Nodes(0,i_aus)
        y5 = Tdomain%Coord_Nodes(1,i_aus)
        z5 = Tdomain%Coord_Nodes(2,i_aus)

        i_aus = Tdomain%specel(n)%Control_Nodes(6)
        x6 = Tdomain%Coord_Nodes(0,i_aus)
        y6 = Tdomain%Coord_Nodes(1,i_aus)
        z6 = Tdomain%Coord_Nodes(2,i_aus)

        i_aus = Tdomain%specel(n)%Control_Nodes(7)
        x7 = Tdomain%Coord_Nodes(0,i_aus)
        y7 = Tdomain%Coord_Nodes(1,i_aus)
        z7 = Tdomain%Coord_Nodes(2,i_aus)

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz


        mat = Tdomain%specel(n)%mat_index

        allocate (Tdomain%specel(n)%Jacob(0:ngllx-1,0:nglly-1,0:ngllz-1) )
        allocate (Tdomain%specel(n)%InvGrad(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2) )


        do k = 0,ngllz - 1
            zeta = Tdomain%sSubdomain(mat)%GLLcz (k)
            do j = 0,nglly - 1
                eta = Tdomain%sSubdomain(mat)%GLLcy (j)
                do i = 0,ngllx - 1
                    xi = Tdomain%sSubdomain(mat)%GLLcx (i)

                    ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)


                    xp = 0.125 * (x0*(1-xi)*(1-eta)*(1-zeta) + x1*(1+xi)*(1-eta)*(1-zeta) + x2*(1+xi)*(1+eta)*(1-zeta) + &
                        x3*(1-xi)*(1+eta)*(1-zeta) + x4*(1-xi)*(1-eta)*(1+zeta) + x5*(1+xi)*(1-eta)*(1+zeta) + &
                        x6*(1+xi)*(1+eta)*(1+zeta) + x7*(1-xi)*(1+eta)*(1+zeta))

                    yp = 0.125 * (y0*(1-xi)*(1-eta)*(1-zeta) + y1*(1+xi)*(1-eta)*(1-zeta) + y2*(1+xi)*(1+eta)*(1-zeta) + &
                        y3*(1-xi)*(1+eta)*(1-zeta) + y4*(1-xi)*(1-eta)*(1+zeta) + y5*(1+xi)*(1-eta)*(1+zeta) + &
                        y6*(1+xi)*(1+eta)*(1+zeta) + y7*(1-xi)*(1+eta)*(1+zeta))

                    zp = 0.125 * (z0*(1-xi)*(1-eta)*(1-zeta) + z1*(1+xi)*(1-eta)*(1-zeta) + z2*(1+xi)*(1+eta)*(1-zeta) + &
                        z3*(1-xi)*(1+eta)*(1-zeta) + z4*(1-xi)*(1-eta)*(1+zeta) + z5*(1+xi)*(1-eta)*(1+zeta) + &
                        z6*(1+xi)*(1+eta)*(1+zeta) + z7*(1-xi)*(1+eta)*(1+zeta))

                    Tdomain%GlobCoord (0,ipoint) = xp
                    Tdomain%GlobCoord (1,ipoint) = yp
                    Tdomain%GlobCoord (2,ipoint) = zp

!!! Computation of the derivative matrix, dx_(jj)/dxi_(ii) !!!
                    LocInvGrad(0,0) = 0.125 * ((x1-x0)*(1-eta)*(1-zeta) + (x2-x3)*(1+eta)*(1-zeta) + &   ! dx/dxi1
                        (x5-x4)*(1-eta)*(1+zeta) + (x6-x7)*(1+eta)*(1+zeta))
                    LocInvGrad(1,0) = 0.125 * ((x3-x0)*(1-xi)*(1-zeta) + (x2-x1)*(1+xi)*(1-zeta) + &     ! dx/dxi2
                        (x7-x4)*(1-xi)*(1+zeta) + (x6-x5)*(1+xi)*(1+zeta))
                    LocInvGrad(2,0) = 0.125 * ((x4-x0)*(1-xi)*(1-eta) + (x5-x1)*(1+xi)*(1-eta) + &       ! dx/dxi3
                        (x7-x3)*(1-xi)*(1+eta) + (x6-x2)*(1+xi)*(1+eta))
                    LocInvGrad(0,1) = 0.125 * ((y1-y0)*(1-eta)*(1-zeta) + (y2-y3)*(1+eta)*(1-zeta) + &   ! dy/dxi1
                        (y5-y4)*(1-eta)*(1+zeta) + (y6-y7)*(1+eta)*(1+zeta))
                    LocInvGrad(1,1) = 0.125 * ((y3-y0)*(1-xi)*(1-zeta) + (y2-y1)*(1+xi)*(1-zeta) + &     ! dy/dxi2
                        (y7-y4)*(1-xi)*(1+zeta) + (y6-y5)*(1+xi)*(1+zeta))
                    LocInvGrad(2,1) = 0.125 * ((y4-y0)*(1-xi)*(1-eta) + (y5-y1)*(1+xi)*(1-eta) + &       ! dy/dxi3
                        (y7-y3)*(1-xi)*(1+eta) + (y6-y2)*(1+xi)*(1+eta))
                    LocInvGrad(0,2) = 0.125 * ((z1-z0)*(1-eta)*(1-zeta) + (z2-z3)*(1+eta)*(1-zeta) + &   ! dz/dxi1
                        (z5-z4)*(1-eta)*(1+zeta) + (z6-z7)*(1+eta)*(1+zeta))
                    LocInvGrad(1,2) = 0.125 * ((z3-z0)*(1-xi)*(1-zeta) + (z2-z1)*(1+xi)*(1-zeta) + &     ! dz/dxi2
                        (z7-z4)*(1-xi)*(1+zeta) + (z6-z5)*(1+xi)*(1+zeta))
                    LocInvGrad(2,2) = 0.125 * ((z4-z0)*(1-xi)*(1-eta) + (z5-z1)*(1+xi)*(1-eta) + &       ! dz/dxi3
                        (z7-z3)*(1-xi)*(1+eta) + (z6-z2)*(1+xi)*(1+eta))

                    call invert_3d (LocInvGrad, Jac)

                    Tdomain%specel(n)%Jacob(i,j,k) = Jac
                    Tdomain%specel(n)%InvGrad (i,j,k, 0:2,0:2) = LocInvGrad(0:2,0:2)
                    !if (abs(Jac) < 0.1 ) write(60,*) rg,n,Jac
                enddo
            enddo
        enddo

    enddo



    ! Super Object Plane Wave
    if (Tdomain%logicD%super_object_local_present) then
        if (Tdomain%super_object_type == "P") then ! Plane Wave

            call semname_shape8_norm(rg,fnamef)
            open (60, file=fnamef, status="unknown", form="formatted")

            call semname_shape8_coord(rg,fnamef)
            open (62, file=fnamef, status="unknown", form="formatted")

            do nf = 0, Tdomain%sPlaneW%n_faces-1

                mat = Tdomain%sPlaneW%pFace(nf)%mat_index

                ngll1 = Tdomain%sPlaneW%pFace(nf)%ngll1
                ngll2 = Tdomain%sPlaneW%pFace(nf)%ngll2

                allocate (Tdomain%sPlaneW%pFace(nf)%normal(0:ngll1-1,0:ngll2-1,0:2))

                if ( Tdomain%sPlaneW%pFace(nf)%dir == 0) then
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(0)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x0 = Tdomain%Coord_Nodes(0,i_aus); y0 = Tdomain%Coord_Nodes(1,i_aus); z0 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(1)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x1 = Tdomain%Coord_Nodes(0,i_aus); y1 = Tdomain%Coord_Nodes(1,i_aus); z1 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(2)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x2 = Tdomain%Coord_Nodes(0,i_aus); y2 = Tdomain%Coord_Nodes(1,i_aus); z2 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(3)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x3 = Tdomain%Coord_Nodes(0,i_aus); y3 = Tdomain%Coord_Nodes(1,i_aus); z3 = Tdomain%Coord_Nodes(2,i_aus)

                    do j = 0,ngll2 - 1
                        eta = Tdomain%sSubdomain(mat)%GLLcy(j)
                        do i = 0,ngll1 - 1
                            xi = Tdomain%sSubdomain(mat)%GLLcx(i)
                            zeta = -1.

                            LocInvGradS(0,0) = 0.125 * ((x1-x0)*(1-eta)*(1-zeta) + (x2-x3)*(1+eta)*(1-zeta) + & ! dx/dxi1
                                (x5-x4)*(1-eta)*(1+zeta) + (x6-x7)*(1+eta)*(1+zeta))
                            LocInvGradS(0,1) = 0.125 * ((x3-x0)*(1-xi)*(1-zeta) + (x2-x1)*(1+xi)*(1-zeta) + &   ! dx/dxi2
                                (x7-x4)*(1-xi)*(1+zeta) + (x6-x5)*(1+xi)*(1+zeta))
                            LocInvGradS(1,0) = 0.125 * ((y1-y0)*(1-eta)*(1-zeta) + (y2-y3)*(1+eta)*(1-zeta) + & ! dy/dxi1
                                (y5-y4)*(1-eta)*(1+zeta) + (y6-y7)*(1+eta)*(1+zeta))
                            LocInvGradS(1,1) = 0.125 * ((y3-y0)*(1-xi)*(1-zeta) + (y2-y1)*(1+xi)*(1-zeta) + &   ! dy/dxi2
                                (y7-y4)*(1-xi)*(1+zeta) + (y6-y5)*(1+xi)*(1+zeta))
                            LocInvGradS(2,0) = 0.125 * ((z1-z0)*(1-eta)*(1-zeta) + (z2-z3)*(1+eta)*(1-zeta) + & ! dz/dxi1
                                (z5-z4)*(1-eta)*(1+zeta) + (z6-z7)*(1+eta)*(1+zeta))
                            LocInvGradS(2,1) = 0.125 * ((z3-z0)*(1-xi)*(1-zeta) + (z2-z1)*(1+xi)*(1-zeta) + &   ! dz/dxi2
                                (z7-z4)*(1-xi)*(1+zeta) + (z6-z5)*(1+xi)*(1+zeta))

                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = -LocInvGradS(0,0)*LocInvGradS(2,1)+LocInvGradS(2,0)*LocInvGradS(0,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                            if ( Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(i,j,0) < 0 ) then
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,0)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,1)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,2)
                            endif

                            if ( Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face0",Tdomain%sPlaneW%pFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(i,j,0)

                        enddo
                    enddo

                else if ( Tdomain%sPlaneW%pFace(nf)%dir == 1) then
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(0)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x0 = Tdomain%Coord_Nodes(0,i_aus); y0 = Tdomain%Coord_Nodes(1,i_aus); z0 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(1)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x1 = Tdomain%Coord_Nodes(0,i_aus); y1 = Tdomain%Coord_Nodes(1,i_aus); z1 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(2)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x5 = Tdomain%Coord_Nodes(0,i_aus); y5 = Tdomain%Coord_Nodes(1,i_aus); z5 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(3)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x4 = Tdomain%Coord_Nodes(0,i_aus); y4 = Tdomain%Coord_Nodes(1,i_aus); z4 = Tdomain%Coord_Nodes(2,i_aus)

                    do j = 0,ngll2 - 1
                        zeta = Tdomain%sSubdomain(mat)%GLLcz(j)
                        do i = 0,ngll1 - 1
                            xi = Tdomain%sSubdomain(mat)%GLLcx(i)
                            eta = -1.

                            LocInvGradS(0,0) = 0.125 * ((x1-x0)*(1-eta)*(1-zeta) + (x2-x3)*(1+eta)*(1-zeta) + & ! dx/dxi1
                                (x5-x4)*(1-eta)*(1+zeta) + (x6-x7)*(1+eta)*(1+zeta))
                            LocInvGradS(0,1) = 0.125 * ((x4-x0)*(1-xi)*(1-eta) + (x5-x1)*(1+xi)*(1-eta) + &     ! dx/dxi3
                                (x7-x3)*(1-xi)*(1+eta) + (x6-x2)*(1+xi)*(1+eta))
                            LocInvGradS(1,0) = 0.125 * ((y1-y0)*(1-eta)*(1-zeta) + (y2-y3)*(1+eta)*(1-zeta) + & ! dy/dxi1
                                (y5-y4)*(1-eta)*(1+zeta) + (y6-y7)*(1+eta)*(1+zeta))
                            LocInvGradS(1,1) = 0.125 * ((y4-y0)*(1-xi)*(1-eta) + (y5-y1)*(1+xi)*(1-eta) + &     ! dy/dxi3
                                (y7-y3)*(1-xi)*(1+eta) + (y6-y2)*(1+xi)*(1+eta))
                            LocInvGradS(2,0) = 0.125 * ((z1-z0)*(1-eta)*(1-zeta) + (z2-z3)*(1+eta)*(1-zeta) + & ! dz/dxi1
                                (z5-z4)*(1-eta)*(1+zeta) + (z6-z7)*(1+eta)*(1+zeta))
                            LocInvGradS(2,1) = 0.125 * ((z4-z0)*(1-xi)*(1-eta) + (z5-z1)*(1+xi)*(1-eta) + &     ! dz/dxi3
                                (z7-z3)*(1-xi)*(1+eta) + (z6-z2)*(1+xi)*(1+eta))

                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = -LocInvGradS(0,0)*LocInvGradS(2,1)+LocInvGradS(2,0)*LocInvGradS(0,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                            if ( Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(i,0,j) > 0 ) then
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,0)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,1)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,2)
                            endif

                            if ( Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face1",Tdomain%sPlaneW%pFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(i,0,j)
                        enddo
                    enddo

                else if ( Tdomain%sPlaneW%pFace(nf)%dir == 2) then
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(0)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x1 = Tdomain%Coord_Nodes(0,i_aus); y1 = Tdomain%Coord_Nodes(1,i_aus); z1 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(1)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x2 = Tdomain%Coord_Nodes(0,i_aus); y2 = Tdomain%Coord_Nodes(1,i_aus); z2 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(2)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x6 = Tdomain%Coord_Nodes(0,i_aus); y6 = Tdomain%Coord_Nodes(1,i_aus); z6 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(3)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x5 = Tdomain%Coord_Nodes(0,i_aus); y5 = Tdomain%Coord_Nodes(1,i_aus); z5 = Tdomain%Coord_Nodes(2,i_aus)

                    do j = 0,ngll2 - 1
                        zeta = Tdomain%sSubdomain(mat)%GLLcz (j)
                        do i = 0,ngll1 - 1
                            eta = Tdomain%sSubdomain(mat)%GLLcy (i)
                            xi = 1.

                            LocInvGradS(0,0) = 0.125 * ((x3-x0)*(1-xi)*(1-zeta) + (x2-x1)*(1+xi)*(1-zeta) + &
                                (x7-x4)*(1-xi)*(1+zeta) + (x6-x5)*(1+xi)*(1+zeta))
                            LocInvGradS(0,1) = 0.125 * ((x4-x0)*(1-xi)*(1-eta) + (x5-x1)*(1+xi)*(1-eta) + &
                                (x7-x3)*(1-xi)*(1+eta) + (x6-x2)*(1+xi)*(1+eta))
                            LocInvGradS(1,0) = 0.125 * ((y3-y0)*(1-xi)*(1-zeta) + (y2-y1)*(1+xi)*(1-zeta) + &
                                (y7-y4)*(1-xi)*(1+zeta) + (y6-y5)*(1+xi)*(1+zeta))
                            LocInvGradS(1,1) = 0.125 * ((y4-y0)*(1-xi)*(1-eta) + (y5-y1)*(1+xi)*(1-eta) + &
                                (y7-y3)*(1-xi)*(1+eta) + (y6-y2)*(1+xi)*(1+eta))
                            LocInvGradS(2,0) = 0.125 * ((z3-z0)*(1-xi)*(1-zeta) + (z2-z1)*(1+xi)*(1-zeta) + &
                                (z7-z4)*(1-xi)*(1+zeta) + (z6-z5)*(1+xi)*(1+zeta))
                            LocInvGradS(2,1) = 0.125 * ((z4-z0)*(1-xi)*(1-eta) + (z5-z1)*(1+xi)*(1-eta) + &
                                (z7-z3)*(1-xi)*(1+eta) + (z6-z2)*(1+xi)*(1+eta))

                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = -LocInvGradS(0,0)*LocInvGradS(2,1)+LocInvGradS(2,0)*LocInvGradS(0,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                            ngll = Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%ngllx

                            if ( Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(ngll-1,i,j) > 0 ) then
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,0)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,1)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,2)
                            endif

                            if ( Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face2",Tdomain%sPlaneW%pFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(5,i,j)

                        enddo
                    enddo

                else if ( Tdomain%sPlaneW%pFace(nf)%dir == 3) then

                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(0)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x3 = Tdomain%Coord_Nodes(0,i_aus); y3 = Tdomain%Coord_Nodes(1,i_aus); z3 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(1)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x2 = Tdomain%Coord_Nodes(0,i_aus); y2 = Tdomain%Coord_Nodes(1,i_aus); z2 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(2)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x6 = Tdomain%Coord_Nodes(0,i_aus); y6 = Tdomain%Coord_Nodes(1,i_aus); z6 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(3)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x7 = Tdomain%Coord_Nodes(0,i_aus); y7 = Tdomain%Coord_Nodes(1,i_aus); z7 = Tdomain%Coord_Nodes(2,i_aus)

                    do j = 0,ngll2 - 1
                        zeta = Tdomain%sSubdomain(mat)%GLLcz (j)
                        do i = 0,ngll1 - 1
                            xi = Tdomain%sSubdomain(mat)%GLLcx (i)
                            eta = 1.

                            LocInvGradS(0,0) = 0.125 * ((x1-x0)*(1-eta)*(1-zeta) + (x2-x3)*(1+eta)*(1-zeta) + & ! dx/dxi1
                                (x5-x4)*(1-eta)*(1+zeta) + (x6-x7)*(1+eta)*(1+zeta))
                            LocInvGradS(0,1) = 0.125 * ((x4-x0)*(1-xi)*(1-eta) + (x5-x1)*(1+xi)*(1-eta) + &     ! dx/dxi3
                                (x7-x3)*(1-xi)*(1+eta) + (x6-x2)*(1+xi)*(1+eta))
                            LocInvGradS(1,0) = 0.125 * ((y1-y0)*(1-eta)*(1-zeta) + (y2-y3)*(1+eta)*(1-zeta) + & ! dy/dxi1
                                (y5-y4)*(1-eta)*(1+zeta) + (y6-y7)*(1+eta)*(1+zeta))
                            LocInvGradS(1,1) = 0.125 * ((y4-y0)*(1-xi)*(1-eta) + (y5-y1)*(1+xi)*(1-eta) + &     ! dy/dxi3
                                (y7-y3)*(1-xi)*(1+eta) + (y6-y2)*(1+xi)*(1+eta))
                            LocInvGradS(2,0) = 0.125 * ((z1-z0)*(1-eta)*(1-zeta) + (z2-z3)*(1+eta)*(1-zeta) + & ! dz/dxi1
                                (z5-z4)*(1-eta)*(1+zeta) + (z6-z7)*(1+eta)*(1+zeta))
                            LocInvGradS(2,1) = 0.125 * ((z4-z0)*(1-xi)*(1-eta) + (z5-z1)*(1+xi)*(1-eta) + &     ! dz/dxi3
                                (z7-z3)*(1-xi)*(1+eta) + (z6-z2)*(1+xi)*(1+eta))

                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = -LocInvGradS(0,0)*LocInvGradS(2,1)+LocInvGradS(2,0)*LocInvGradS(0,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                            ngll = Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%nglly

                            if ( Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(i,ngll-1,j) < 0 ) then
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,0)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,1)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,2)
                            endif

                            if ( Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face3",Tdomain%sPlaneW%pFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(i,5,j)
                        enddo
                    enddo

                else if ( Tdomain%sPlaneW%pFace(nf)%dir == 4) then

                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(0)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x0 = Tdomain%Coord_Nodes(0,i_aus); y0 = Tdomain%Coord_Nodes(1,i_aus); z0 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(1)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x3 = Tdomain%Coord_Nodes(0,i_aus); y3 = Tdomain%Coord_Nodes(1,i_aus); z3 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(2)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x7 = Tdomain%Coord_Nodes(0,i_aus); y7 = Tdomain%Coord_Nodes(1,i_aus); z7 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(3)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x4 = Tdomain%Coord_Nodes(0,i_aus); y4 = Tdomain%Coord_Nodes(1,i_aus); z4 = Tdomain%Coord_Nodes(2,i_aus)

                    do j = 0,ngll2 - 1
                        zeta = Tdomain%sSubdomain(mat)%GLLcz (j)
                        do i = 0,ngll1 - 1
                            eta = Tdomain%sSubdomain(mat)%GLLcy (i)
                            xi = -1.

                            LocInvGradS(0,0) = 0.125 * ((x3-x0)*(1-xi)*(1-zeta) + (x2-x1)*(1+xi)*(1-zeta) + &
                                (x7-x4)*(1-xi)*(1+zeta) + (x6-x5)*(1+xi)*(1+zeta))
                            LocInvGradS(0,1) = 0.125 * ((x4-x0)*(1-xi)*(1-eta) + (x5-x1)*(1+xi)*(1-eta) + &
                                (x7-x3)*(1-xi)*(1+eta) + (x6-x2)*(1+xi)*(1+eta))
                            LocInvGradS(1,0) = 0.125 * ((y3-y0)*(1-xi)*(1-zeta) + (y2-y1)*(1+xi)*(1-zeta) + &
                                (y7-y4)*(1-xi)*(1+zeta) + (y6-y5)*(1+xi)*(1+zeta))
                            LocInvGradS(1,1) = 0.125 * ((y4-y0)*(1-xi)*(1-eta) + (y5-y1)*(1+xi)*(1-eta) + &
                                (y7-y3)*(1-xi)*(1+eta) + (y6-y2)*(1+xi)*(1+eta))
                            LocInvGradS(2,0) = 0.125 * ((z3-z0)*(1-xi)*(1-zeta) + (z2-z1)*(1+xi)*(1-zeta) + &
                                (z7-z4)*(1-xi)*(1+zeta) + (z6-z5)*(1+xi)*(1+zeta))
                            LocInvGradS(2,1) = 0.125 * ((z4-z0)*(1-xi)*(1-eta) + (z5-z1)*(1+xi)*(1-eta) + &
                                (z7-z3)*(1-xi)*(1+eta) + (z6-z2)*(1+xi)*(1+eta))

                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = -LocInvGradS(0,0)*LocInvGradS(2,1)+LocInvGradS(2,0)*LocInvGradS(0,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                            if ( Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(0,i,j) < 0 ) then
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,0)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,1)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,2)
                            endif

                            if ( Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face4",Tdomain%sPlaneW%pFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(0,i,j)

                        enddo
                    enddo

                else if ( Tdomain%sPlaneW%pFace(nf)%dir == 5) then

                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(0)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x4 = Tdomain%Coord_Nodes(0,i_aus); y4 = Tdomain%Coord_Nodes(1,i_aus); z4 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(1)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x5 = Tdomain%Coord_Nodes(0,i_aus); y5 = Tdomain%Coord_Nodes(1,i_aus); z5 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(2)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x6 = Tdomain%Coord_Nodes(0,i_aus); y6 = Tdomain%Coord_Nodes(1,i_aus); z6 = Tdomain%Coord_Nodes(2,i_aus)
                    nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(3)
                    i_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP; x7 = Tdomain%Coord_Nodes(0,i_aus); y7 = Tdomain%Coord_Nodes(1,i_aus); z7 = Tdomain%Coord_Nodes(2,i_aus)

                    do j = 0,ngll2 - 1
                        eta = Tdomain%sSubdomain(mat)%GLLcy (j)
                        do i = 0,ngll1 - 1
                            xi = Tdomain%sSubdomain(mat)%GLLcx (i)
                            zeta = 1.

                            LocInvGradS(0,0) = 0.125 * ((x1-x0)*(1-eta)*(1-zeta) + (x2-x3)*(1+eta)*(1-zeta) + & ! dx/dxi1
                                (x5-x4)*(1-eta)*(1+zeta) + (x6-x7)*(1+eta)*(1+zeta))
                            LocInvGradS(0,1) = 0.125 * ((x3-x0)*(1-xi)*(1-zeta) + (x2-x1)*(1+xi)*(1-zeta) + &   ! dx/dxi2
                                (x7-x4)*(1-xi)*(1+zeta) + (x6-x5)*(1+xi)*(1+zeta))
                            LocInvGradS(1,0) = 0.125 * ((y1-y0)*(1-eta)*(1-zeta) + (y2-y3)*(1+eta)*(1-zeta) + & ! dy/dxi1
                                (y5-y4)*(1-eta)*(1+zeta) + (y6-y7)*(1+eta)*(1+zeta))
                            LocInvGradS(1,1) = 0.125 * ((y3-y0)*(1-xi)*(1-zeta) + (y2-y1)*(1+xi)*(1-zeta) + &   ! dy/dxi2
                                (y7-y4)*(1-xi)*(1+zeta) + (y6-y5)*(1+xi)*(1+zeta))
                            LocInvGradS(2,0) = 0.125 * ((z1-z0)*(1-eta)*(1-zeta) + (z2-z3)*(1+eta)*(1-zeta) + & ! dz/dxi1
                                (z5-z4)*(1-eta)*(1+zeta) + (z6-z7)*(1+eta)*(1+zeta))
                            LocInvGradS(2,1) = 0.125 * ((z3-z0)*(1-xi)*(1-zeta) + (z2-z1)*(1+xi)*(1-zeta) + &   ! dz/dxi2
                                (z7-z4)*(1-xi)*(1+zeta) + (z6-z5)*(1+xi)*(1+zeta))

                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = -LocInvGradS(0,0)*LocInvGradS(2,1)+LocInvGradS(2,0)*LocInvGradS(0,1)
                            Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                            ngll = Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%ngllz

                            if ( Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(i,j,ngll-1) > 0 ) then
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,0) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,0)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,1) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,1)
                                Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) = - Tdomain%sPlaneW%pFace(nf)%normal(i,j,2)
                            endif

                            if ( Tdomain%sPlaneW%pFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face5",Tdomain%sPlaneW%pFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sPlaneW%pFace(nf)%Face_UP)%Which_Elem)%Jacob(i,j,5)

                        enddo
                    enddo
                endif

                allocate (Tdomain%sPlaneW%pFace(nf)%Coord_nodes(1:ngll1-2,1:ngll2-2,0:2))

                do j = 1,ngll2-2
                    do i = 1,ngll1-2
                        n_aus = Tdomain%sPlaneW%pFace(nf)%Face_UP
                        Tdomain%sPlaneW%pFace(nf)%Coord_nodes(i,j,0:2) = Tdomain%GlobCoord(0:2,Tdomain%sFace(n_aus)%Iglobnum_Face(i,j))
                    enddo
                enddo

                !!
                !! Test pour verifier si l'orientation entre les faces UP et DOWN est bonne - base sur les coordonnees !!
                !
                n_up = Tdomain%sPlaneW%pFace(nf)%Face_UP
                n_down = Tdomain%sPlaneW%pFace(nf)%Face_DOWN
                if ( Tdomain%sPlaneW%pFace(nf)%Orient == 0 ) then
                    do j = 1,ngll2-2
                        do i = 1,ngll1-2
                            if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(0,Tdomain%sFace(n_down)%Iglobnum_Face(i,j)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(1,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(1,Tdomain%sFace(n_down)%Iglobnum_Face(i,j)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(2,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(2,Tdomain%sFace(n_down)%Iglobnum_Face(i,j)) )<0.1  ) ) then
                                write(62,*) 'Face',nv,Tdomain%sPlaneW%pEdge(ne)%Orient
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sPlaneW%pFace(nf)%Orient == 1 ) then
                    do j = 1,ngll2-2
                        do i = 1,ngll1-2
                            if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(0,Tdomain%sFace(n_down)%Iglobnum_Face(ngll1-1-i,j)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(1,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(1,Tdomain%sFace(n_down)%Iglobnum_Face(ngll1-1-i,j)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(2,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(2,Tdomain%sFace(n_down)%Iglobnum_Face(ngll1-1-i,j)) )<0.1  ) ) then
                                write(62,*) 'Face',nv,Tdomain%sPlaneW%pEdge(ne)%Orient
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sPlaneW%pFace(nf)%Orient == 2 ) then
                    do j = 1,ngll2-2
                        do i = 1,ngll1-2
                            if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(0,Tdomain%sFace(n_down)%Iglobnum_Face(i,ngll2-1-j)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(1,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(1,Tdomain%sFace(n_down)%Iglobnum_Face(i,ngll2-1-j)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(2,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(2,Tdomain%sFace(n_down)%Iglobnum_Face(i,ngll2-1-j)) )<0.1  ) ) then
                                write(62,*) 'Face',nv,Tdomain%sPlaneW%pEdge(ne)%Orient
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sPlaneW%pFace(nf)%Orient == 3 ) then
                    do j = 1,ngll2-2
                        do i = 1,ngll1-2
                            if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(0,Tdomain%sFace(n_down)%Iglobnum_Face(ngll1-1-i,ngll2-1-j)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(1,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(1,Tdomain%sFace(n_down)%Iglobnum_Face(ngll1-1-i,ngll2-1-j)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(2,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(2,Tdomain%sFace(n_down)%Iglobnum_Face(ngll1-1-i,ngll2-1-j)) )<0.1  ) ) then
                                write(62,*) 'Face',nv,Tdomain%sPlaneW%pEdge(ne)%Orient
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sPlaneW%pFace(nf)%Orient == 4 ) then
                    do j = 1,ngll2-2
                        do i = 1,ngll1-2
                            if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(0,Tdomain%sFace(n_down)%Iglobnum_Face(j,i)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(1,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(1,Tdomain%sFace(n_down)%Iglobnum_Face(j,i)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(2,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(2,Tdomain%sFace(n_down)%Iglobnum_Face(j,i)) )<0.1  ) ) then
                                write(62,*) 'Face',nv,Tdomain%sPlaneW%pEdge(ne)%Orient
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sPlaneW%pFace(nf)%Orient == 5 ) then
                    do j = 1,ngll2-2
                        do i = 1,ngll1-2
                            if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(0,Tdomain%sFace(n_down)%Iglobnum_Face(ngll2-1-j,i)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(1,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(1,Tdomain%sFace(n_down)%Iglobnum_Face(ngll2-1-j,i)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(2,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(2,Tdomain%sFace(n_down)%Iglobnum_Face(ngll2-1-j,i)) )<0.1  ) ) then
                                write(62,*) 'Face',nv,Tdomain%sPlaneW%pEdge(ne)%Orient
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sPlaneW%pFace(nf)%Orient == 6 ) then
                    do j = 1,ngll2-2
                        do i = 1,ngll1-2
                            if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(0,Tdomain%sFace(n_down)%Iglobnum_Face(j,ngll1-1-i)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(1,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(1,Tdomain%sFace(n_down)%Iglobnum_Face(j,ngll1-1-i)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(2,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(2,Tdomain%sFace(n_down)%Iglobnum_Face(j,ngll1-1-i)) )<0.1  ) ) then
                                write(62,*) 'Face',nv,Tdomain%sPlaneW%pEdge(ne)%Orient
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sPlaneW%pFace(nf)%Orient == 7 ) then
                    do j = 1,ngll2-2
                        do i = 1,ngll1-2
                            if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(0,Tdomain%sFace(n_down)%Iglobnum_Face(ngll2-1-j,ngll1-1-i)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(1,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(1,Tdomain%sFace(n_down)%Iglobnum_Face(ngll2-1-j,ngll1-1-i)) )<0.1  .and. &
                                abs( Tdomain%GlobCoord(2,Tdomain%sFace(n_up)%Iglobnum_Face(i,j)) - Tdomain%GlobCoord(2,Tdomain%sFace(n_down)%Iglobnum_Face(ngll2-1-j,ngll1-1-i)) )<0.1  ) ) then
                                write(62,*) 'Face',nv,Tdomain%sPlaneW%pEdge(ne)%Orient
                            endif
                        enddo
                    enddo
                endif
                !!
            enddo

            do ne = 0, Tdomain%sPlaneW%n_edges-1
                ngll = Tdomain%sPlaneW%pEdge(ne)%ngll
                allocate (Tdomain%sPlaneW%pEdge(ne)%Coord_nodes(1:ngll-2,0:2))
                do i = 1,ngll-2
                    n_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_UP
                    Tdomain%sPlaneW%pEdge(ne)%Coord_nodes(i,0:2) = Tdomain%GlobCoord(0:2,Tdomain%sEdge(n_aus)%Iglobnum_Edge(i))
                enddo
                !
                !! Test pour verifier si l'orientation entre les edges UP et DOWN est bonne - base sur les coordonnees !!
                !
                n_up = Tdomain%sPlaneW%pEdge(ne)%Edge_UP
                n_down = Tdomain%sPlaneW%pEdge(ne)%Edge_DOWN
                if ( Tdomain%sPlaneW%pEdge(ne)%Orient == 0 ) then
                    do i = 1,ngll-2
                        if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sEdge(n_up)%Iglobnum_Edge(i)) - Tdomain%GlobCoord(0,Tdomain%sEdge(n_down)%Iglobnum_Edge(i)) )<0.1  .and. &
                            abs( Tdomain%GlobCoord(1,Tdomain%sEdge(n_up)%Iglobnum_Edge(i)) - Tdomain%GlobCoord(1,Tdomain%sEdge(n_down)%Iglobnum_Edge(i)) )<0.1  .and. &
                            abs( Tdomain%GlobCoord(2,Tdomain%sEdge(n_up)%Iglobnum_Edge(i)) - Tdomain%GlobCoord(2,Tdomain%sEdge(n_down)%Iglobnum_Edge(i)) )<0.1  ) ) then
                            write(62,*) 'Edge',nv,Tdomain%sPlaneW%pEdge(ne)%Orient
                        endif
                    enddo
                else
                    do i = 1,ngll-2
                        if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sEdge(n_up)%Iglobnum_Edge(i)) - Tdomain%GlobCoord(0,Tdomain%sEdge(n_down)%Iglobnum_Edge(ngll-1-i)) )<0.1  .and. &
                            abs( Tdomain%GlobCoord(1,Tdomain%sEdge(n_up)%Iglobnum_Edge(i)) - Tdomain%GlobCoord(1,Tdomain%sEdge(n_down)%Iglobnum_Edge(ngll-1-i)) )<0.1  .and. &
                            abs( Tdomain%GlobCoord(2,Tdomain%sEdge(n_up)%Iglobnum_Edge(i)) - Tdomain%GlobCoord(2,Tdomain%sEdge(n_down)%Iglobnum_Edge(ngll-1-i)) )<0.1  ) ) then
                            write(62,*) 'Edge',nv,Tdomain%sPlaneW%pEdge(ne)%Orient
                        endif
                    enddo
                endif
                !!
            enddo

            do nv = 0, Tdomain%sPlaneW%n_vertices-1
                n_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP
                Tdomain%sPlaneW%pVertex(nv)%Coord_nodes(0:2) = Tdomain%GlobCoord(0:2,Tdomain%sVertex(n_aus)%Iglobnum_Vertex)

                !
                !! Test pour verifier si l'orientation entre les vertex UP et DOWN est bonne - base sur les coordonnees !!
                !
                n_up = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP
                n_down = Tdomain%sPlaneW%pVertex(nv)%Vertex_DOWN
                if (.not. ( abs( Tdomain%GlobCoord(0,Tdomain%sVertex(n_up)%Iglobnum_Vertex) - Tdomain%GlobCoord(0,Tdomain%sVertex(n_down)%Iglobnum_Vertex) )<0.1  .and. &
                    abs( Tdomain%GlobCoord(1,Tdomain%sVertex(n_up)%Iglobnum_Vertex) - Tdomain%GlobCoord(1,Tdomain%sVertex(n_down)%Iglobnum_Vertex) )<0.1  .and. &
                    abs( Tdomain%GlobCoord(2,Tdomain%sVertex(n_up)%Iglobnum_Vertex) - Tdomain%GlobCoord(2,Tdomain%sVertex(n_down)%Iglobnum_Vertex) )<0.1  ) ) then
                    write(62,*) 'Vertex',nv
                endif
                !!
            enddo

            close(60)
            close(62)

        endif
    endif

    !stop


    ! Neumann
    if (Tdomain%logicD%neumann_local_present) then

        call semname_shape8_neu(rg,fnamef)
        open (60, file=fnamef, status="unknown", form="formatted")

        do nf = 0, Tdomain%sNeu%n_faces-1


            mat = Tdomain%sNeu%nFace(nf)%mat_index
            ngll1 = Tdomain%sNeu%nFace(nf)%ngll1
            ngll2 = Tdomain%sNeu%nFace(nf)%ngll2

            !    allocate (Tdomain%sNeu%nFace(nf)%ds(0:ngll1-1,0:ngll2-1))
            allocate (Tdomain%sNeu%nFace(nf)%normal(0:ngll1-1,0:ngll2-1,0:2))

            if ( Tdomain%sNeu%nFace(nf)%dir == 0) then
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(0)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x0 = Tdomain%Coord_Nodes(0,i_aus); y0 = Tdomain%Coord_Nodes(1,i_aus); z0 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(1)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x1 = Tdomain%Coord_Nodes(0,i_aus); y1 = Tdomain%Coord_Nodes(1,i_aus); z1 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(2)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x2 = Tdomain%Coord_Nodes(0,i_aus); y2 = Tdomain%Coord_Nodes(1,i_aus); z2 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(3)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x3 = Tdomain%Coord_Nodes(0,i_aus); y3 = Tdomain%Coord_Nodes(1,i_aus); z3 = Tdomain%Coord_Nodes(2,i_aus)


                do j = 0,ngll2 - 1
                    eta = Tdomain%sSubdomain(mat)%GLLcy(j)
                    do i = 0,ngll1 - 1
                        xi = Tdomain%sSubdomain(mat)%GLLcx(i)
                        zeta = -1.

                        LocInvGradS(0,0) = 0.125 * ((x1-x0)*(1-eta)*(1-zeta) + (x2-x3)*(1+eta)*(1-zeta) + & ! dx/dxi1
                            (x5-x4)*(1-eta)*(1+zeta) + (x6-x7)*(1+eta)*(1+zeta))
                        LocInvGradS(0,1) = 0.125 * ((x3-x0)*(1-xi)*(1-zeta) + (x2-x1)*(1+xi)*(1-zeta) + & ! dx/dxi2
                            (x7-x4)*(1-xi)*(1+zeta) + (x6-x5)*(1+xi)*(1+zeta))
                        LocInvGradS(1,0) = 0.125 * ((y1-y0)*(1-eta)*(1-zeta) + (y2-y3)*(1+eta)*(1-zeta) + & ! dy/dxi1
                            (y5-y4)*(1-eta)*(1+zeta) + (y6-y7)*(1+eta)*(1+zeta))
                        LocInvGradS(1,1) = 0.125 * ((y3-y0)*(1-xi)*(1-zeta) + (y2-y1)*(1+xi)*(1-zeta) + & ! dy/dxi2
                            (y7-y4)*(1-xi)*(1+zeta) + (y6-y5)*(1+xi)*(1+zeta))
                        LocInvGradS(2,0) = 0.125 * ((z1-z0)*(1-eta)*(1-zeta) + (z2-z3)*(1+eta)*(1-zeta) + & ! dz/dxi1
                            (z5-z4)*(1-eta)*(1+zeta) + (z6-z7)*(1+eta)*(1+zeta))
                        LocInvGradS(2,1) = 0.125 * ((z3-z0)*(1-xi)*(1-zeta) + (z2-z1)*(1+xi)*(1-zeta) + & ! dz/dxi2
                            (z7-z4)*(1-xi)*(1+zeta) + (z6-z5)*(1+xi)*(1+zeta))

                        Tdomain%sNeu%nFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,1) = LocInvGradS(2,0)*LocInvGradS(0,1)-LocInvGradS(0,0)*LocInvGradS(2,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                        norm = (LocInvGradS(0,0)**2+LocInvGradS(1,0)**2+LocInvGradS(2,0)**2)*(LocInvGradS(0,1)**2+LocInvGradS(1,1)**2+LocInvGradS(2,1)**2)- &
                            (LocInvGradS(0,0)*LocInvGradS(0,1)+LocInvGradS(1,0)*LocInvGradS(1,1)+LocInvGradS(2,0)*LocInvGradS(2,1))**2

                        if ( Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(i,j,0) > 0 ) then
                            Tdomain%sNeu%nFace(nf)%normal(i,j,0) = - Tdomain%sNeu%nFace(nf)%normal(i,j,0)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,1) = - Tdomain%sNeu%nFace(nf)%normal(i,j,1)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,2) = - Tdomain%sNeu%nFace(nf)%normal(i,j,2)
                        endif

                        if ( Tdomain%sNeu%nFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face0",Tdomain%sNeu%nFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(i,j,0)

                    enddo
                enddo

            else if ( Tdomain%sNeu%nFace(nf)%dir == 1) then

                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(0)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x0 = Tdomain%Coord_Nodes(0,i_aus); y0 = Tdomain%Coord_Nodes(1,i_aus); z0 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(1)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x1 = Tdomain%Coord_Nodes(0,i_aus); y1 = Tdomain%Coord_Nodes(1,i_aus); z1 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(2)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x5 = Tdomain%Coord_Nodes(0,i_aus); y5 = Tdomain%Coord_Nodes(1,i_aus); z5 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(3)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x4 = Tdomain%Coord_Nodes(0,i_aus); y4 = Tdomain%Coord_Nodes(1,i_aus); z4 = Tdomain%Coord_Nodes(2,i_aus)

                do j = 0,ngll2 - 1
                    zeta = Tdomain%sSubdomain(mat)%GLLcz(j)
                    do i = 0,ngll1 - 1
                        xi = Tdomain%sSubdomain(mat)%GLLcx(i)
                        eta = -1.

                        LocInvGradS(0,0) = 0.125 * ((x1-x0)*(1-eta)*(1-zeta) + (x2-x3)*(1+eta)*(1-zeta) + & ! dx/dxi1
                            (x5-x4)*(1-eta)*(1+zeta) + (x6-x7)*(1+eta)*(1+zeta))
                        LocInvGradS(0,1) = 0.125 * ((x4-x0)*(1-xi)*(1-eta) + (x5-x1)*(1+xi)*(1-eta) + &     ! dx/dxi3
                            (x7-x3)*(1-xi)*(1+eta) + (x6-x2)*(1+xi)*(1+eta))
                        LocInvGradS(1,0) = 0.125 * ((y1-y0)*(1-eta)*(1-zeta) + (y2-y3)*(1+eta)*(1-zeta) + & ! dy/dxi1
                            (y5-y4)*(1-eta)*(1+zeta) + (y6-y7)*(1+eta)*(1+zeta))
                        LocInvGradS(1,1) = 0.125 * ((y4-y0)*(1-xi)*(1-eta) + (y5-y1)*(1+xi)*(1-eta) + &     ! dy/dxi3
                            (y7-y3)*(1-xi)*(1+eta) + (y6-y2)*(1+xi)*(1+eta))
                        LocInvGradS(2,0) = 0.125 * ((z1-z0)*(1-eta)*(1-zeta) + (z2-z3)*(1+eta)*(1-zeta) + & ! dz/dxi1
                            (z5-z4)*(1-eta)*(1+zeta) + (z6-z7)*(1+eta)*(1+zeta))
                        LocInvGradS(2,1) = 0.125 * ((z4-z0)*(1-xi)*(1-eta) + (z5-z1)*(1+xi)*(1-eta) + &     ! dz/dxi3
                            (z7-z3)*(1-xi)*(1+eta) + (z6-z2)*(1+xi)*(1+eta))

                        Tdomain%sNeu%nFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,1) = LocInvGradS(2,0)*LocInvGradS(0,1)-LocInvGradS(0,0)*LocInvGradS(2,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                        norm = (LocInvGradS(0,0)**2+LocInvGradS(1,0)**2+LocInvGradS(2,0)**2)*(LocInvGradS(0,1)**2+LocInvGradS(1,1)**2+LocInvGradS(2,1)**2)- &
                            (LocInvGradS(0,0)*LocInvGradS(0,1)+LocInvGradS(1,0)*LocInvGradS(1,1)+LocInvGradS(2,0)*LocInvGradS(2,1))**2

                        if ( Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(i,0,j) < 0 ) then
                            Tdomain%sNeu%nFace(nf)%normal(i,j,0) = - Tdomain%sNeu%nFace(nf)%normal(i,j,0)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,1) = - Tdomain%sNeu%nFace(nf)%normal(i,j,1)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,2) = - Tdomain%sNeu%nFace(nf)%normal(i,j,2)
                        endif

                        if ( Tdomain%sNeu%nFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face1",Tdomain%sNeu%nFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(i,0,j)

                    enddo
                enddo

            else if ( Tdomain%sNeu%nFace(nf)%dir == 2) then

                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(0)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x1 = Tdomain%Coord_Nodes(0,i_aus); y1 = Tdomain%Coord_Nodes(1,i_aus); z1 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(1)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x2 = Tdomain%Coord_Nodes(0,i_aus); y2 = Tdomain%Coord_Nodes(1,i_aus); z2 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(2)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x6 = Tdomain%Coord_Nodes(0,i_aus); y6 = Tdomain%Coord_Nodes(1,i_aus); z6 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(3)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x5 = Tdomain%Coord_Nodes(0,i_aus); y5 = Tdomain%Coord_Nodes(1,i_aus); z5 = Tdomain%Coord_Nodes(2,i_aus)

                do j = 0,ngll2 - 1
                    zeta = Tdomain%sSubdomain(mat)%GLLcz (j)
                    do i = 0,ngll1 - 1
                        eta = Tdomain%sSubdomain(mat)%GLLcy (i)
                        xi = 1.

                        LocInvGradS(0,0) = 0.125 * ((x3-x0)*(1-xi)*(1-zeta) + (x2-x1)*(1+xi)*(1-zeta) + &
                            (x7-x4)*(1-xi)*(1+zeta) + (x6-x5)*(1+xi)*(1+zeta))
                        LocInvGradS(0,1) = 0.125 * ((x4-x0)*(1-xi)*(1-eta) + (x5-x1)*(1+xi)*(1-eta) + &
                            (x7-x3)*(1-xi)*(1+eta) + (x6-x2)*(1+xi)*(1+eta))
                        LocInvGradS(1,0) = 0.125 * ((y3-y0)*(1-xi)*(1-zeta) + (y2-y1)*(1+xi)*(1-zeta) + &
                            (y7-y4)*(1-xi)*(1+zeta) + (y6-y5)*(1+xi)*(1+zeta))
                        LocInvGradS(1,1) = 0.125 * ((y4-y0)*(1-xi)*(1-eta) + (y5-y1)*(1+xi)*(1-eta) + &
                            (y7-y3)*(1-xi)*(1+eta) + (y6-y2)*(1+xi)*(1+eta))
                        LocInvGradS(2,0) = 0.125 * ((z3-z0)*(1-xi)*(1-zeta) + (z2-z1)*(1+xi)*(1-zeta) + &
                            (z7-z4)*(1-xi)*(1+zeta) + (z6-z5)*(1+xi)*(1+zeta))
                        LocInvGradS(2,1) = 0.125 * ((z4-z0)*(1-xi)*(1-eta) + (z5-z1)*(1+xi)*(1-eta) + &
                            (z7-z3)*(1-xi)*(1+eta) + (z6-z2)*(1+xi)*(1+eta))

                        Tdomain%sNeu%nFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,1) = LocInvGradS(2,0)*LocInvGradS(0,1)-LocInvGradS(0,0)*LocInvGradS(2,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                        norm = (LocInvGradS(0,0)**2+LocInvGradS(1,0)**2+LocInvGradS(2,0)**2)*(LocInvGradS(0,1)**2+LocInvGradS(1,1)**2+LocInvGradS(2,1)**2)- &
                            (LocInvGradS(0,0)*LocInvGradS(0,1)+LocInvGradS(1,0)*LocInvGradS(1,1)+LocInvGradS(2,0)*LocInvGradS(2,1))**2

                        ngll = Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%ngllx

                        if ( Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(ngll-1,i,j) < 0 ) then
                            Tdomain%sNeu%nFace(nf)%normal(i,j,0) = - Tdomain%sNeu%nFace(nf)%normal(i,j,0)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,1) = - Tdomain%sNeu%nFace(nf)%normal(i,j,1)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,2) = - Tdomain%sNeu%nFace(nf)%normal(i,j,2)
                        endif

                        if ( Tdomain%sNeu%nFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face2",nf,Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem,Tdomain%sNeu%nFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(5,i,j)

                    enddo
                enddo

            else if ( Tdomain%sNeu%nFace(nf)%dir == 3) then

                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(0)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x3 = Tdomain%Coord_Nodes(0,i_aus); y3 = Tdomain%Coord_Nodes(1,i_aus); z3 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(1)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x2 = Tdomain%Coord_Nodes(0,i_aus); y2 = Tdomain%Coord_Nodes(1,i_aus); z2 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(2)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x6 = Tdomain%Coord_Nodes(0,i_aus); y6 = Tdomain%Coord_Nodes(1,i_aus); z6 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(3)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x7 = Tdomain%Coord_Nodes(0,i_aus); y7 = Tdomain%Coord_Nodes(1,i_aus); z7 = Tdomain%Coord_Nodes(2,i_aus)

                do j = 0,ngll2 - 1
                    zeta = Tdomain%sSubdomain(mat)%GLLcz (j)
                    do i = 0,ngll1 - 1
                        xi = Tdomain%sSubdomain(mat)%GLLcx (i)
                        eta = 1.

                        LocInvGradS(0,0) = 0.125 * ((x1-x0)*(1-eta)*(1-zeta) + (x2-x3)*(1+eta)*(1-zeta) + & ! dx/dxi1
                            (x5-x4)*(1-eta)*(1+zeta) + (x6-x7)*(1+eta)*(1+zeta))
                        LocInvGradS(0,1) = 0.125 * ((x4-x0)*(1-xi)*(1-eta) + (x5-x1)*(1+xi)*(1-eta) + &     ! dx/dxi3
                            (x7-x3)*(1-xi)*(1+eta) + (x6-x2)*(1+xi)*(1+eta))
                        LocInvGradS(1,0) = 0.125 * ((y1-y0)*(1-eta)*(1-zeta) + (y2-y3)*(1+eta)*(1-zeta) + & ! dy/dxi1
                            (y5-y4)*(1-eta)*(1+zeta) + (y6-y7)*(1+eta)*(1+zeta))
                        LocInvGradS(1,1) = 0.125 * ((y4-y0)*(1-xi)*(1-eta) + (y5-y1)*(1+xi)*(1-eta) + &     ! dy/dxi3
                            (y7-y3)*(1-xi)*(1+eta) + (y6-y2)*(1+xi)*(1+eta))
                        LocInvGradS(2,0) = 0.125 * ((z1-z0)*(1-eta)*(1-zeta) + (z2-z3)*(1+eta)*(1-zeta) + & ! dz/dxi1
                            (z5-z4)*(1-eta)*(1+zeta) + (z6-z7)*(1+eta)*(1+zeta))
                        LocInvGradS(2,1) = 0.125 * ((z4-z0)*(1-xi)*(1-eta) + (z5-z1)*(1+xi)*(1-eta) + &     ! dz/dxi3
                            (z7-z3)*(1-xi)*(1+eta) + (z6-z2)*(1+xi)*(1+eta))

                        Tdomain%sNeu%nFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,1) = LocInvGradS(2,0)*LocInvGradS(0,1)-LocInvGradS(0,0)*LocInvGradS(2,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                        norm = (LocInvGradS(0,0)**2+LocInvGradS(1,0)**2+LocInvGradS(2,0)**2)*(LocInvGradS(0,1)**2+LocInvGradS(1,1)**2+LocInvGradS(2,1)**2)- &
                            (LocInvGradS(0,0)*LocInvGradS(0,1)+LocInvGradS(1,0)*LocInvGradS(1,1)+LocInvGradS(2,0)*LocInvGradS(2,1))**2

                        ngll = Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%nglly

                        if ( Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(i,ngll-1,j) > 0 ) then
                            Tdomain%sNeu%nFace(nf)%normal(i,j,0) = - Tdomain%sNeu%nFace(nf)%normal(i,j,0)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,1) = - Tdomain%sNeu%nFace(nf)%normal(i,j,1)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,2) = - Tdomain%sNeu%nFace(nf)%normal(i,j,2)
                        endif

                        if ( Tdomain%sNeu%nFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face3",Tdomain%sNeu%nFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(i,5,j)

                    enddo
                enddo

            else if ( Tdomain%sNeu%nFace(nf)%dir == 4) then

                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(0)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x0 = Tdomain%Coord_Nodes(0,i_aus); y0 = Tdomain%Coord_Nodes(1,i_aus); z0 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(1)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x3 = Tdomain%Coord_Nodes(0,i_aus); y3 = Tdomain%Coord_Nodes(1,i_aus); z3 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(2)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x7 = Tdomain%Coord_Nodes(0,i_aus); y7 = Tdomain%Coord_Nodes(1,i_aus); z7 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(3)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x4 = Tdomain%Coord_Nodes(0,i_aus); y4 = Tdomain%Coord_Nodes(1,i_aus); z4 = Tdomain%Coord_Nodes(2,i_aus)

                do j = 0,ngll2 - 1
                    zeta = Tdomain%sSubdomain(mat)%GLLcz (j)
                    do i = 0,ngll1 - 1
                        eta = Tdomain%sSubdomain(mat)%GLLcy (i)
                        xi = -1.

                        LocInvGradS(0,0) = 0.125 * ((x3-x0)*(1-xi)*(1-zeta) + (x2-x1)*(1+xi)*(1-zeta) + &
                            (x7-x4)*(1-xi)*(1+zeta) + (x6-x5)*(1+xi)*(1+zeta))
                        LocInvGradS(0,1) = 0.125 * ((x4-x0)*(1-xi)*(1-eta) + (x5-x1)*(1+xi)*(1-eta) + &
                            (x7-x3)*(1-xi)*(1+eta) + (x6-x2)*(1+xi)*(1+eta))
                        LocInvGradS(1,0) = 0.125 * ((y3-y0)*(1-xi)*(1-zeta) + (y2-y1)*(1+xi)*(1-zeta) + &
                            (y7-y4)*(1-xi)*(1+zeta) + (y6-y5)*(1+xi)*(1+zeta))
                        LocInvGradS(1,1) = 0.125 * ((y4-y0)*(1-xi)*(1-eta) + (y5-y1)*(1+xi)*(1-eta) + &
                            (y7-y3)*(1-xi)*(1+eta) + (y6-y2)*(1+xi)*(1+eta))
                        LocInvGradS(2,0) = 0.125 * ((z3-z0)*(1-xi)*(1-zeta) + (z2-z1)*(1+xi)*(1-zeta) + &
                            (z7-z4)*(1-xi)*(1+zeta) + (z6-z5)*(1+xi)*(1+zeta))
                        LocInvGradS(2,1) = 0.125 * ((z4-z0)*(1-xi)*(1-eta) + (z5-z1)*(1+xi)*(1-eta) + &
                            (z7-z3)*(1-xi)*(1+eta) + (z6-z2)*(1+xi)*(1+eta))

                        Tdomain%sNeu%nFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,1) = LocInvGradS(2,0)*LocInvGradS(0,1)-LocInvGradS(0,0)*LocInvGradS(2,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                        norm = (LocInvGradS(0,0)**2+LocInvGradS(1,0)**2+LocInvGradS(2,0)**2)*(LocInvGradS(0,1)**2+LocInvGradS(1,1)**2+LocInvGradS(2,1)**2)- &
                            (LocInvGradS(0,0)*LocInvGradS(0,1)+LocInvGradS(1,0)*LocInvGradS(1,1)+LocInvGradS(2,0)*LocInvGradS(2,1))**2

                        if ( Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(0,i,j) > 0 ) then
                            Tdomain%sNeu%nFace(nf)%normal(i,j,0) = - Tdomain%sNeu%nFace(nf)%normal(i,j,0)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,1) = - Tdomain%sNeu%nFace(nf)%normal(i,j,1)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,2) = - Tdomain%sNeu%nFace(nf)%normal(i,j,2)
                        endif

                        if ( Tdomain%sNeu%nFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face4",Tdomain%sNeu%nFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(0,i,j)

                    enddo
                enddo

            else if ( Tdomain%sNeu%nFace(nf)%dir == 5) then

                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(0)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x4 = Tdomain%Coord_Nodes(0,i_aus); y4 = Tdomain%Coord_Nodes(1,i_aus); z4 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(1)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x5 = Tdomain%Coord_Nodes(0,i_aus); y5 = Tdomain%Coord_Nodes(1,i_aus); z5 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(2)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x6 = Tdomain%Coord_Nodes(0,i_aus); y6 = Tdomain%Coord_Nodes(1,i_aus); z6 = Tdomain%Coord_Nodes(2,i_aus)
                nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(3)
                i_aus = Tdomain%sNeu%nVertex(nv)%Vertex; x7 = Tdomain%Coord_Nodes(0,i_aus); y7 = Tdomain%Coord_Nodes(1,i_aus); z7 = Tdomain%Coord_Nodes(2,i_aus)

                do j = 0,ngll2 - 1
                    eta = Tdomain%sSubdomain(mat)%GLLcy (j)
                    do i = 0,ngll1 - 1
                        xi = Tdomain%sSubdomain(mat)%GLLcx (i)
                        zeta = 1.

                        LocInvGradS(0,0) = 0.125 * ((x1-x0)*(1-eta)*(1-zeta) + (x2-x3)*(1+eta)*(1-zeta) + & ! dx/dxi1
                            (x5-x4)*(1-eta)*(1+zeta) + (x6-x7)*(1+eta)*(1+zeta))
                        LocInvGradS(0,1) = 0.125 * ((x3-x0)*(1-xi)*(1-zeta) + (x2-x1)*(1+xi)*(1-zeta) + & ! dx/dxi2
                            (x7-x4)*(1-xi)*(1+zeta) + (x6-x5)*(1+xi)*(1+zeta))
                        LocInvGradS(1,0) = 0.125 * ((y1-y0)*(1-eta)*(1-zeta) + (y2-y3)*(1+eta)*(1-zeta) + & ! dy/dxi1
                            (y5-y4)*(1-eta)*(1+zeta) + (y6-y7)*(1+eta)*(1+zeta))
                        LocInvGradS(1,1) = 0.125 * ((y3-y0)*(1-xi)*(1-zeta) + (y2-y1)*(1+xi)*(1-zeta) + & ! dy/dxi2
                            (y7-y4)*(1-xi)*(1+zeta) + (y6-y5)*(1+xi)*(1+zeta))
                        LocInvGradS(2,0) = 0.125 * ((z1-z0)*(1-eta)*(1-zeta) + (z2-z3)*(1+eta)*(1-zeta) + & ! dz/dxi1
                            (z5-z4)*(1-eta)*(1+zeta) + (z6-z7)*(1+eta)*(1+zeta))
                        LocInvGradS(2,1) = 0.125 * ((z3-z0)*(1-xi)*(1-zeta) + (z2-z1)*(1+xi)*(1-zeta) + & ! dz/dxi2
                            (z7-z4)*(1-xi)*(1+zeta) + (z6-z5)*(1+xi)*(1+zeta))

                        Tdomain%sNeu%nFace(nf)%normal(i,j,0) = LocInvGradS(1,0)*LocInvGradS(2,1)-LocInvGradS(2,0)*LocInvGradS(1,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,1) = LocInvGradS(2,0)*LocInvGradS(0,1)-LocInvGradS(0,0)*LocInvGradS(2,1)
                        Tdomain%sNeu%nFace(nf)%normal(i,j,2) = LocInvGradS(0,0)*LocInvGradS(1,1)-LocInvGradS(1,0)*LocInvGradS(0,1)

                        norm = ( LocInvGradS(0,0)**2 + LocInvGradS(1,0)**2 + LocInvGradS(2,0)**2 )*( LocInvGradS(0,1)**2 + LocInvGradS(1,1)**2 + LocInvGradS(2,1)**2 ) -&
                            ( LocInvGradS(0,0)*LocInvGradS(0,1) + LocInvGradS(1,0)*LocInvGradS(1,1) + LocInvGradS(2,0)*LocInvGradS(2,1) )**2

                        ngll = Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%ngllz

                        if ( Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(i,j,ngll-1) < 0 ) then
                            Tdomain%sNeu%nFace(nf)%normal(i,j,0) = - Tdomain%sNeu%nFace(nf)%normal(i,j,0)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,1) = - Tdomain%sNeu%nFace(nf)%normal(i,j,1)
                            Tdomain%sNeu%nFace(nf)%normal(i,j,2) = - Tdomain%sNeu%nFace(nf)%normal(i,j,2)
                        endif

                        if ( Tdomain%sNeu%nFace(nf)%normal(i,j,2) < 0 ) write(60,*) "Face5",Tdomain%sNeu%nFace(nf)%normal(i,j,2),Tdomain%specel(Tdomain%sFace(Tdomain%sNeu%nFace(nf)%Face)%Which_Elem)%Jacob(i,j,5)

                    enddo
                enddo

            endif

            allocate (Tdomain%sNeu%nFace(nf)%Coord_nodes(1:ngll1-2,1:ngll2-2,0:2))
            do j = 1,ngll2-2
                do i = 1,ngll1-2
                    n_aus = Tdomain%sNeu%nFace(nf)%Face
                    Tdomain%sNeu%nFace(nf)%Coord_nodes(i,j,0:2) = Tdomain%GlobCoord(0:2,Tdomain%sFace(n_aus)%Iglobnum_Face(i,j))
                enddo
            enddo

        enddo

        do ne = 0, Tdomain%sNeu%n_edges-1
            ngll = Tdomain%sNeu%nEdge(ne)%ngll
            allocate (Tdomain%sNeu%nEdge(ne)%Coord_nodes(1:ngll-2,0:2))
            do i = 1,ngll-2
                n_aus = Tdomain%sNeu%nEdge(ne)%Edge
                Tdomain%sNeu%nEdge(ne)%Coord_nodes(i,0:2) = Tdomain%GlobCoord(0:2,Tdomain%sEdge(n_aus)%Iglobnum_Edge(i))
            enddo
        enddo

        do nv = 0, Tdomain%sNeu%n_vertices-1
            n_aus = Tdomain%sNeu%nVertex(nv)%Vertex
            Tdomain%sNeu%nVertex(nv)%Coord_nodes(0:2) = Tdomain%GlobCoord(0:2,Tdomain%sVertex(n_aus)%Iglobnum_Vertex)
        enddo

        close(60)

    endif ! Neumann


    do n = 0,Tdomain%n_elem - 1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        do k = 0,ngllz - 1
            do j = 0,nglly - 1
                do i = 0,ngllx - 1
                    Tdomain%specel(n)%Jacob(i,j,k) = abs(Tdomain%specel(n)%Jacob(i,j,k))
                enddo
            enddo
        enddo
    enddo

    !!purge - fuites memoire
    do n = 0,Tdomain%n_face-1
        deallocate (Tdomain%sFace(n)%Iglobnum_Face)
    enddo

    do n = 0,Tdomain%n_edge-1
        deallocate (Tdomain%sEdge(n)%Iglobnum_Edge)
    enddo

    return
end subroutine shape8
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
