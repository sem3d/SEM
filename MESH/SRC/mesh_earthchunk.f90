module mesh_earthchunk

    implicit none

    real, parameter :: Earth_radius=6371000.d0, Pi = 3.141592653, Pi180 = Pi/180.0
    type :: EarthChunk_t
        integer :: npt_lon, npt_lat, npt_z, total_pt, total_elem, nods
        real :: lon_center, lat_center, delta_lon, delta_lat,&
                  max_depth, step_lonlat
        real,  dimension(:), allocatable :: z_depth
    end type



contains

    subroutine sph2cart(r, lon, lat, x, y, z)
        implicit none
        real, intent(in) :: r, lon, lat
        real, intent(out) :: x, y, z
        real :: colatRad, lonRad
        colatRad = (90.0-lat)*Pi180
        lonRad = lon*Pi180

        x = r * sin(colatRad) * sin(lonRad)
        y = r * sin(colatRad) * cos(lonRad)
        z = r * cos(colatRad)


    end subroutine sph2cart


    function find_nextdepth(curDepth, delta_angle, ratio) result(z)
        implicit none
        real :: z
        real, intent(in) :: delta_angle, curDepth, ratio

        z = delta_angle * (Earth_radius-curDepth) * ratio + curDepth
    end function find_nextdepth


    subroutine init_earthchunk(chunk)
        implicit none
        real :: delta_lon, delta_lat, step_lonlat, max_depth, curDepth, dTheta
        integer :: mesh_type, countPt, i
        type(EarthChunk_t), intent(out) :: chunk

        write(*,*) "*****************************************"
        write(*,*) "  --> Longitude center: " ; read(*,*) chunk%lon_center
        write(*,*) "  --> Latitude center : " ; read(*,*) chunk%lat_center
        write(*,*) "  --> Total length longitude: " ; read(*,*) delta_lon
        write(*,*) "  --> Total length latitude : " ; read(*,*) delta_lat
        write(*,*) "  --> Spatial step (in deg) : " ; read(*,*) step_lonlat
        write(*,*) "  --> Max depth (in meter)  : " ; read(*,*) chunk%max_depth
        write(*,*) "  --> 8 or 27 nodes (1 or 2)?: " ; read(*,*) mesh_type


!        if(mesh_type /= 1 .and. mesh_type /= 2) &
!            stop "init_earthchunk : mesh_type = [12]"
!        if(xmin > xmax .or. ymin > ymax .or. zmin > zmax .or. step_x <= 0d0 .or. step_y <= 0d0 .or. step_z <= 0d0) &
!            stop "init_earthchunk : bad limits"


        chunk%step_lonlat = step_lonlat

        chunk%npt_lon = ceiling(delta_lon/chunk%step_lonlat)+1+2 ! +2 pour les pml
        chunk%npt_lat = ceiling(delta_lat/chunk%step_lonlat)+1+2 ! +2 pour les pml
        chunk%delta_lon = (chunk%npt_lon-1) * chunk%step_lonlat
        chunk%delta_lat = (chunk%npt_lat-1) * chunk%step_lonlat



        countPt=0
        dTheta = Pi/180.0 *chunk%step_lonlat
        curDepth=0.0
        do while(curDepth < chunk%max_depth)
            curDepth = find_nextdepth(curDepth, dTheta, 0.5)
            countPt=countPt+1
        end do
        chunk%npt_z = countPt+1 ! +1 pour la pml du fond



        chunk%total_elem = (chunk%npt_lon-1) * (chunk%npt_lat-1) * (chunk%npt_z-1)



        if( mesh_type == 2) then
            chunk%nods=27
            chunk%npt_lon = 2*chunk%npt_lon -1
            chunk%npt_lat = 2*chunk%npt_lat -1
            chunk%npt_z = 2*chunk%npt_z-1

            allocate(chunk%z_depth(0:chunk%npt_z))

            curDepth=0.0
            chunk%z_depth(0) = curDepth
            do i=1,chunk%npt_z-1,2
                curDepth = find_nextdepth(curDepth, dTheta, 0.5)
                chunk%z_depth(i)   = 0.5*(chunk%z_depth(i-1) + curDepth)
                chunk%z_depth(i+1) = curDepth
            end do

            chunk%step_lonlat = step_lonlat/2.0


        else
            chunk%nods=8

            allocate(chunk%z_depth(0:chunk%npt_z))

            curDepth=0.0
            chunk%z_depth(0) = curDepth
            do i=1,chunk%npt_z-1
                curDepth = find_nextdepth(curDepth, dTheta, 0.5)
                chunk%z_depth(i) = curDepth
            end do

        endif

        chunk%total_pt = chunk%npt_lon * chunk%npt_lat * chunk%npt_z


!        write(*,*) "pt ",chunk%total_pt
!        write(*,*) "elem ",chunk%total_elem
!        write(*,*) "delta lon ",chunk%delta_lon, "delta lat ", chunk%delta_lat 

    end subroutine init_earthchunk


    subroutine clean_earthchunk(chunk)
        implicit none
        type(EarthChunk_t), intent(inout) :: chunk

        deallocate(chunk%z_depth)




    end subroutine clean_earthchunk



    subroutine defineRotation(lon, lat, mat)
        real, intent(in) :: lon, lat
        real, dimension(0:2,0:2) :: mat
        real :: theta, phi

        theta = lat*Pi180
        phi = lon*Pi180
        mat(0,0) = cos(phi)
        mat(0,1) = sin(phi)*cos(theta)
        mat(0,2) = -sin(phi)*sin(theta)
        mat(1,0) = -sin(phi)
        mat(1,1) = cos(phi)*cos(theta)
        mat(1,2) = -cos(phi)*sin(theta)
        mat(2,0) = 0
        mat(2,1) = sin(theta)
        mat(2,2) = cos(theta)

    end subroutine defineRotation


    subroutine create_earthchunk(chunk, pml_b, nmat, xp, yp, zp, Ipoint, mat)
        implicit none
        integer, intent(in) :: nmat, pml_b
        real, dimension(0:), intent(out) :: xp, yp, zp
        integer, dimension(0:,0:), intent(out) :: Ipoint
        integer, dimension(0:), intent(out) :: mat
        integer :: nelemx, nelemy, nelemz, iex, iey, iez, ix, iy, iz, indice, dX, dXY
        real :: curDepth, curX, curY, x, y, z
        real, dimension(0:2,0:2) :: rotToRealChunk
!        real :: refX, refY, sizeChunk
        type(EarthChunk_t), intent(in) :: chunk


        call defineRotation(chunk%lon_center, chunk%lat_center, rotToRealChunk)



        indice=0
        do iz=0,chunk%npt_z-1
            curDepth = chunk%z_depth(iz)
            do iy=0, chunk%npt_lon-1
                curY = - chunk%delta_lon/2 + iy*chunk%step_lonlat
                do ix=0, chunk%npt_lat-1
                    curX = - chunk%delta_lat/2 + ix*chunk%step_lonlat
                    call sph2cart(Earth_radius-curDepth, curY, curX, x, y, z)
                    xp(indice) = x * rotToRealChunk(0,0) + y * rotToRealChunk(0,1) + z * rotToRealChunk(0,2)
                    yp(indice) = x * rotToRealChunk(1,0) + y * rotToRealChunk(1,1) + z * rotToRealChunk(1,2)
                    zp(indice) = x * rotToRealChunk(2,0) + y * rotToRealChunk(2,1) + z * rotToRealChunk(2,2)
                    indice = indice+1
                end do
            end do
        end do

        indice=0
            
        nelemx = chunk%npt_lat-1
        nelemy = chunk%npt_lon-1
        nelemz = chunk%npt_z-1
        dX = chunk%npt_lat
        dXY = chunk%npt_lat*chunk%npt_lon

        if( chunk%nods == 27 ) then
            nelemx = (chunk%npt_lat-1)/2
            nelemy = (chunk%npt_lon-1)/2
            nelemz = (chunk%npt_z-1)/2
        endif

        do iez=0,nelemz-1
            do iey=0, nelemy-1
                do iex=0, nelemx-1

                    if( chunk%nods == 27 ) then
                        ix = 2*iex
                        iy = 2*iey
                        iz = 2*iez

                        Ipoint( 0, indice) = ix   + iy    *dX + iz    *dXY
                        Ipoint( 1, indice) = ix+2 + iy    *dX + iz    *dXY
                        Ipoint( 2, indice) = ix+2 + (iy+2)*dX + iz    *dXY
                        Ipoint( 3, indice) = ix   + (iy+2)*dX + iz    *dXY
                        Ipoint( 4, indice) = ix   + iy    *dX + (iz+2)*dXY
                        Ipoint( 5, indice) = ix+2 + iy    *dX + (iz+2)*dXY
                        Ipoint( 6, indice) = ix+2 + (iy+2)*dX + (iz+2)*dXY
                        Ipoint( 7, indice) = ix   + (iy+2)*dX + (iz+2)*dXY

                        Ipoint( 8, indice) = ix+1 + iy    *dX + iz    *dXY
                        Ipoint( 9, indice) = ix+2 + (iy+1)*dX + iz    *dXY
                        Ipoint(10, indice) = ix+1 + (iy+2)*dX + iz    *dXY
                        Ipoint(11, indice) = ix   + (iy+1)*dX + iz    *dXY

                        Ipoint(12, indice) = ix   + iy    *dX + (iz+1)*dXY
                        Ipoint(13, indice) = ix+2 + iy    *dX + (iz+1)*dXY
                        Ipoint(14, indice) = ix+2 + (iy+2)*dX + (iz+1)*dXY
                        Ipoint(15, indice) = ix   + (iy+2)*dX + (iz+1)*dXY

                        Ipoint(16, indice) = ix+1 + iy    *dX + (iz+2)*dXY
                        Ipoint(17, indice) = ix+2 + (iy+1)*dX + (iz+2)*dXY
                        Ipoint(18, indice) = ix+1 + (iy+2)*dX + (iz+2)*dXY
                        Ipoint(19, indice) = ix   + (iy+1)*dX + (iz+2)*dXY

                        Ipoint(20, indice) = ix+1 + (iy+1)*dX + iz    *dXY
                        Ipoint(21, indice) = ix+1 + iy    *dX + (iz+1)*dXY
                        Ipoint(22, indice) = ix+2 + (iy+1)*dX + (iz+1)*dXY
                        Ipoint(23, indice) = ix+1 + (iy+2)*dX + (iz+1)*dXY
                        Ipoint(24, indice) = ix   + (iy+1)*dX + (iz+1)*dXY
                        Ipoint(25, indice) = ix+1 + (iy+1)*dX + (iz+2)*dXY

                        Ipoint(26, indice) = ix+1 + (iy+1)*dX + (iz+1)*dXY

                    else
                        ix = iex
                        iy = iey
                        iz = iez
                        Ipoint(0, indice) = ix   + iy    *dX + iz    *dXY
                        Ipoint(1, indice) = ix+1 + iy    *dX + iz    *dXY
                        Ipoint(2, indice) = ix+1 + (iy+1)*dX + iz    *dXY
                        Ipoint(3, indice) = ix   + (iy+1)*dX + iz    *dXY
                        Ipoint(4, indice) = ix   + iy    *dX + (iz+1)*dXY
                        Ipoint(5, indice) = ix+1 + iy    *dX + (iz+1)*dXY
                        Ipoint(6, indice) = ix+1 + (iy+1)*dX + (iz+1)*dXY
                        Ipoint(7, indice) = ix   + (iy+1)*dX + (iz+1)*dXY

                    endif

                    mat(indice) = nmat-1

                    if( iez == nelemz-1) then
                        mat(indice) = nmat+16
                        if( iex==0 ) then
                            mat(indice) = nmat+15
                            if( iey==0) then
                                mat(indice) = nmat+8
                            else if( iey==nelemy-1) then
                                mat(indice) = nmat+11
                            end if
                        else if( iex==nelemx-1) then
                            mat(indice) = nmat+13
                            if( iey==0) then
                                mat(indice) = nmat+9
                            else if( iey==nelemy-1) then
                                mat(indice) = nmat+10
                            end if
                        else if( iey==0) then
                            mat(indice) = nmat+12
                        else if( iey==nelemy-1) then
                            mat(indice) = nmat+14
                        end if
                    else if( iex==0 ) then
                        mat(indice) = nmat+7
                        if( iey==0) then
                            mat(indice) = nmat+0
                        else if( iey==nelemy-1) then
                            mat(indice) = nmat+3
                        end if
                    else if( iex==nelemx-1) then
                        mat(indice) = nmat+5
                        if( iey==0) then
                            mat(indice) = nmat+1
                        else if( iey==nelemy-1) then
                            mat(indice) = nmat+2
                        end if
                    else if( iey==0 ) then
                        mat(indice) = nmat+4
                    else if( iey== nelemy-1) then
                        mat(indice) = nmat+6
                    end if

                    indice = indice+1
                end do
            end do
        end do

    end subroutine create_earthchunk




end module mesh_earthchunk




        

