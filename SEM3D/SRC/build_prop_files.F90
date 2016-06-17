module build_prop_files

    use sdomain
    use mpi
    use constants
    use scomm
    use scommutils
    use define_random
    use interpolation

    implicit none
#include "index.h"

contains
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine create_prop_files(Tdomain, rg)

        !INPUTS
        type (domain), intent (INOUT), target :: Tdomain
        integer      , intent(IN) :: rg

        !LOCAL
        integer :: mat, assocMat, code


        call build_random_properties(Tdomain, rg)


    end subroutine create_prop_files

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine apply_prop_files(Tdomain, rg)

        use interpolation

        implicit none
        !INPUTS
        type (domain), intent (INOUT), target :: Tdomain
        integer      , intent(IN) :: rg

        !LOCAL
        integer :: mat, n, ipoint, i, j, k, lnum, bnum, ee
        integer :: ngll, elem_mat
        double precision, dimension(:, :), allocatable :: interpolatedRF !Properties
        double precision, dimension(:), allocatable :: lambda
        integer :: assocMat, propId

        !Putting properties on elements
        if(Tdomain%any_Random) then
            allocate(interpolatedRF(0:size(Tdomain%GlobCoord, 2) - 1, 0:nProp - 1))
            allocate(lambda(0:size(Tdomain%GlobCoord, 2) - 1))
        end if

        do mat = 0, Tdomain%n_mat - 1
        !do mat = 0, 0 !For Tests

            write(*,*) "  "
            write(*,*) "  Analyzing Material ", mat, "------------------- in rg ", rg
            !write(*,*) "  Tdomain%subD_exist ", Tdomain%subD_exist, "------------------- in rg ", rg

            if(.not. Tdomain%subD_exist(mat)) cycle

            !write(*,*) "  Exists in proc ", rg

            assocMat = Tdomain%sSubdomain(mat)%assocMat

            if(propOnFile(Tdomain, assocMat)) then
                !write(*,*) "  Material ", mat, " have properties on file"
                if(is_rand(Tdomain%sSubdomain(assocMat))) then

                    do propId = 0, nProp - 1
                    !write(*,*) "  Reading and Interpolating Random Properties on:"
                    !do propId = 0, 0 !For Tests
                        !write(*,*) " prop = ", propId, "rg = ", rg
                        write(*,*) trim(Tdomain%sSubDomain(mat)%propFilePath(propId))
                        call interpolateToMesh_V2(BBoxFileName=Tdomain%sSubDomain(mat)%propFilePath(propId),  &
                                               coordList=Tdomain%GlobCoord, &
                                               UNV_randField=interpolatedRF(:,propId:propId), &
                                               rang=rg, &
                                               xMinLoc_In = Tdomain%sSubDomain(mat)%MinBound_Loc,  &
                                               xMaxLoc_In = Tdomain%sSubDomain(mat)%MaxBound_Loc,  &
                                               mat = mat, &
                                               Tdomain  = Tdomain)
                        !write(*,*) " AFTER prop = ", propId, "rg = ", rg
                    end do

                    !write(*,*) " AFTER 2 rg = ", rg

                    !INFO
                    !interpolatedRF(:, 0) --Density
                    !interpolatedRF(:, 1) --Kappa/Lambda
                    !interpolatedRF(:, 2) --Mu

                    if(Tdomain%sSubDomain(mat)%lambdaSwitch == 1) then
                        lambda(:) = interpolatedRF(:, 1)
                        interpolatedRF(:, 1) = lambda(:) + (2.0d0*interpolatedRF(:, 2)/3.0D0)!Kappa
                    else
                        lambda(:) = interpolatedRF(:, 1) - (2.0D0*interpolatedRF(:, 2)/3.0D0)
                    end if
                    !S%DMu = S%Sspeed**2 * S%Ddensity
                    !S%DLambda = (S%Pspeed**2 - 2 * S%Sspeed **2 ) * S%Ddensity
                    !S%DKappa = S%DLambda + 2.*S%DMu /3.
                    !S%DLambda = S%DKappa - 2.*S%DMu /3.

                    !Applying properties to elements
                    do n = 0, Tdomain%n_elem-1
                        elem_mat = Tdomain%specel(n)%mat_index

                        if(elem_mat /= mat) cycle

                        lnum = Tdomain%specel(n)%lnum
                        bnum = lnum/VCHUNK
                        ee = mod(lnum,VCHUNK)
                        ngll = Tdomain%sSubDomain(mat)%NGLL

                        !Properties by element
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID)
                                if (Tdomain%sdom%n_sls>0)  then
                                    if (Tdomain%sdom%aniso) then
                                        Tdomain%sdom%Q_(:,:,:,bnum,ee) = Tdomain%sSubDomain(mat)%Qmu
                                    else
                                        Tdomain%sdom%Qs_(:,:,:,bnum,ee) = Tdomain%sSubDomain(mat)%Qmu
                                        Tdomain%sdom%Qp_(:,:,:,bnum,ee) = Tdomain%sSubDomain(mat)%Qpression
                                    endif
                                endif
                            case (DM_FLUID)
                                !Nothing to do, all definitions are on GLLs
                            case (DM_SOLID_PML)
                                !Nothing to do, all definitions are on GLLs
                            case (DM_FLUID_PML)
                                !Nothing to do, all definitions are on GLLs
                        end select

                        !Properties by GLL
                        do i = 0, ngll-1
                            do j = 0, ngll-1
                                do k = 0, ngll-1
                                    ipoint  = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    select case (Tdomain%specel(n)%domain)
                                        case (DM_SOLID)
                                            Tdomain%sdom%Density_(i,j,k,bnum,ee) = interpolatedRF(ipoint, 0)
                                            Tdomain%sdom%Lambda_ (i,j,k,bnum,ee) = lambda(ipoint)
                                            Tdomain%sdom%Mu_     (i,j,k,bnum,ee) = interpolatedRF(ipoint, 2)
                                            Tdomain%sdom%Kappa_  (i,j,k,bnum,ee) = interpolatedRF(ipoint, 1)
                                        case (DM_SOLID_PML)
                                            Tdomain%spmldom%Density_(i,j,k,bnum,ee) = interpolatedRF(ipoint, 0)
#ifndef CPML
                                            Tdomain%spmldom%Lambda_ (i,j,k,bnum,ee) = lambda(ipoint)
                                            Tdomain%spmldom%Mu_     (i,j,k,bnum,ee) = interpolatedRF(ipoint, 2)
#endif
                                        case (DM_FLUID)
                                            Tdomain%fdom%IDensity_(i,j,k,bnum,ee) = 1d0/interpolatedRF(ipoint, 0)
                                            Tdomain%fdom%Lambda_  (i,j,k,bnum,ee) = lambda(ipoint)
                                        case (DM_FLUID_PML)
                                            Tdomain%fpmldom%Density_(i,j,k,bnum,ee) = interpolatedRF(ipoint, 0)
                                            Tdomain%fpmldom%Lambda_ (i,j,k,bnum,ee) = lambda(ipoint)
                                    end select
                                end do
                            end do
                        end do
                    end do
                    !write(*,*) " AFTER 3 rg = ", rg
                end if
                !write(*,*) " AFTER 4 rg = ", rg
            end if
            !write(*,*) " AFTER 5 rg = ", rg
        end do

        !write(*,*) " AFTER 6 rg = ", rg

        !Adjusting PML values


        if(allocated(interpolatedRF)) deallocate(interpolatedRF)
        if(allocated(lambda))         deallocate(lambda)


        write(*,*) "  Propagating properties to PML in rg ", rg
        call propagate_PML_properties(Tdomain)

        !write(*,*) " AFTER 7 rg = ", rg

    end subroutine apply_prop_files

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    function propOnFile(Tdomain, mat) result(authorization)

        !INPUTS
        type (domain), intent (in), target :: Tdomain
        integer, intent(in) :: mat

        !OUTPUT
        logical :: authorization

        !LOCAL
        integer :: assocMat

        assocMat = Tdomain%sSubdomain(mat)%assocMat
        authorization = .false.

        if(Tdomain%sSubdomain(assocMat)%initial_material_type == 'R') then
            authorization = .true.
        end if
        !end if

    end function

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine propagate_PML_properties(Tdomain)
        !This subroutine modifies the "prop" table.
        !In the PML elements it extrudes the properties according to the PML direction

        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain

        !LOCAL
        integer :: ngll, lnum, bnum, ee
        integer :: iPoint
        integer :: i, j, k, n !counters
        integer :: LimPML1, LimPML2, LimPML3
        integer :: dir, mat_index
        double precision :: pointDens, pointLambda, pointMu;
        logical :: verbose = .false.

        !Propagating Properties over the PML

        ! With the new PML definition, we can do better:
        ! We can project x,y,z, on the pml plane, and interpolate the properties from the
        ! projected coordinates

        write(*,*) "Tdomain%not_PML_List = ", Tdomain%not_PML_List

        do n = 0, Tdomain%n_elem-1
            mat_index = Tdomain%specel(n)%mat_index

            if(Tdomain%not_PML_List(mat_index)) cycle

            ngll = Tdomain%sSubDomain(mat_index)%NGLL
            lnum = Tdomain%specel(n)%lnum
            bnum = lnum/VCHUNK
            ee = mod(lnum,VCHUNK)
            !On the lower bound initialization
            LimPML1 = 0
            LimPML2 = 0
            LimPML3 = 0
            !To the upper bound transformation
            if (Tdomain%sSubDomain(mat_index)%pml_width(0)<0) LimPML1 = ngll-1
            if (Tdomain%sSubDomain(mat_index)%pml_width(1)<0) LimPML2 = ngll-1
            if (Tdomain%sSubDomain(mat_index)%pml_width(2)<0) LimPML3 = ngll-1

            dir = read_PML_Direction(Tdomain, mat_index)
            if(verbose) write(*,*) "dir = ", dir
            if(verbose) write(*,*) "Tdomain%specel(n)%domain = ", Tdomain%specel(n)%domain
            pointDens = 0

            !FOR TESTS
            !select case (Tdomain%specel(n)%domain)
            !    case (DM_SOLID_PML)
            !        Tdomain%spmldom%Density_(:,:,:,bnum,ee) = 0
            !    case (DM_FLUID_PML)
            !
            !end select



            select case(dir)
                !Face X oriented
                case(0)
                    if(verbose) write(*,*) "Face X oriented"
                    do j = 0, ngll-1
                        do k = 0, ngll-1
                            select case (Tdomain%specel(n)%domain)
                                case (DM_SOLID_PML)
                                    pointDens   = Tdomain%spmldom%Density_(LimPML1,j,k,bnum,ee)
#ifndef CPML
                                    pointLambda = Tdomain%spmldom%Lambda_(LimPML1,j,k,bnum,ee)
                                    pointMu     = Tdomain%spmldom%Mu_(LimPML1,j,k,bnum,ee)
#endif
                                case (DM_FLUID_PML)

                            end select

                            do i = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,bnum,ee) = pointDens
#ifndef CPML
                                        Tdomain%spmldom%Lambda_(i,j,k,bnum,ee)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,bnum,ee)      = pointMu
#endif
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        enddo
                    enddo
                !
                !Face Y oriented
                case(1)
                    if(verbose) write(*,*) "Face Y oriented"
                    do i = 0, ngll-1
                        do k = 0, ngll-1
                            select case (Tdomain%specel(n)%domain)
                                case (DM_SOLID_PML)
                                    pointDens   = Tdomain%spmldom%Density_(i, LimPML2,k,bnum,ee)
#ifndef CPML
                                    pointLambda = Tdomain%spmldom%Lambda_(i, LimPML2,k,bnum,ee)
                                    pointMu     = Tdomain%spmldom%Mu_(i, LimPML2,k,bnum,ee)
#endif
                                case (DM_FLUID_PML)

                            end select

                            do j = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,bnum,ee) = pointDens
#ifndef CPML
                                        Tdomain%spmldom%Lambda_(i,j,k,bnum,ee)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,bnum,ee)      = pointMu
#endif
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        enddo
                    enddo

                !Face Z oriented
                case(2)
                    if(verbose) write(*,*) "Face Z oriented"
                    do i = 0, ngll-1
                        do j = 0, ngll-1
                            select case (Tdomain%specel(n)%domain)
                                case (DM_SOLID_PML)
                                    pointDens   = Tdomain%spmldom%Density_(i,j,LimPML3,bnum,ee)
#ifndef CPML
                                    pointLambda = Tdomain%spmldom%Lambda_(i,j,LimPML3,bnum,ee)
                                    pointMu     = Tdomain%spmldom%Mu_(i,j,LimPML3,bnum,ee)
#endif
                                case (DM_FLUID_PML)

                            end select

                            do k = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,bnum,ee) = pointDens
#ifndef CPML
                                        Tdomain%spmldom%Lambda_(i,j,k,bnum,ee)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,bnum,ee)      = pointMu
#endif
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        enddo
                    enddo

                !Edge in XY
                case(3)
                    if(verbose) write(*,*) "Edge in XY"
                    do k = 0, ngll-1
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_PML)
                                pointDens   = Tdomain%spmldom%Density_(LimPML1,LimPML2,k,bnum,ee)
#ifndef CPML
                                pointLambda = Tdomain%spmldom%Lambda_(LimPML1,LimPML2,k,bnum,ee)
                                pointMu     = Tdomain%spmldom%Mu_(LimPML1,LimPML2,k,bnum,ee)
#endif
                            case (DM_FLUID_PML)

                        end select

                        do i = 0, ngll-1
                            do j = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,bnum,ee) = pointDens
#ifndef CPML
                                        Tdomain%spmldom%Lambda_(i,j,k,bnum,ee)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,bnum,ee)      = pointMu
#endif
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        end do
                    enddo

                !Edge in YZ
                case(4)
                    if(verbose) write(*,*) "Edge in YZ"
                    do i = 0, ngll-1
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_PML)
                                pointDens   = Tdomain%spmldom%Density_(i, LimPML2,LimPML3,bnum,ee)
#ifndef CPML
                                pointLambda = Tdomain%spmldom%Lambda_(i, LimPML2,LimPML3,bnum,ee)
                                pointMu     = Tdomain%spmldom%Mu_(i, LimPML2,LimPML3,bnum,ee)
#endif
                            case (DM_FLUID_PML)

                        end select

                        do j = 0, ngll-1
                            do k = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,bnum,ee) = pointDens
#ifndef CPML
                                        Tdomain%spmldom%Lambda_(i,j,k,bnum,ee)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,bnum,ee)      = pointMu
#endif
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        end do
                    enddo

                !Edge in ZX
                case(5)
                    if(verbose) write(*,*) "Edge in ZX"
                    do j = 0, ngll-1
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_PML)
                                pointDens   = Tdomain%spmldom%Density_(LimPML1, j,LimPML3,bnum,ee)
#ifndef CPML
                                pointLambda = Tdomain%spmldom%Lambda_(LimPML1, j,LimPML3,bnum,ee)
                                pointMu     = Tdomain%spmldom%Mu_(LimPML1, j,LimPML3,bnum,ee)
#endif
                            case (DM_FLUID_PML)

                        end select

                        do i = 0, ngll-1
                            do k = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,bnum,ee) = pointDens
#ifndef CPML
                                        Tdomain%spmldom%Lambda_(i,j,k,bnum,ee)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,bnum,ee)      = pointMu
#endif
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        end do
                    enddo

                !Vertex
                case(6)
                    if(verbose) write(*,*) "Vertex"
                    select case (Tdomain%specel(n)%domain)
                        case (DM_SOLID_PML)
                            pointDens   = Tdomain%spmldom%Density_(LimPML1, LimPML2,LimPML3,bnum,ee)
#ifndef CPML
                            pointLambda = Tdomain%spmldom%Lambda_(LimPML1, LimPML2,LimPML3,bnum,ee)
                            pointMu     = Tdomain%spmldom%Mu_(LimPML1, LimPML2,LimPML3,bnum,ee)
#endif
                        case (DM_FLUID_PML)

                    end select

                    do i = 0, ngll-1
                        do j = 0, ngll-1
                            do k = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,bnum,ee) = pointDens
#ifndef CPML
                                        Tdomain%spmldom%Lambda_(i,j,k,bnum,ee)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,bnum,ee)      = pointMu
#endif
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        end do
                    end do
            end select
        end do

    end subroutine propagate_PML_properties

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    function read_PML_Direction(Tdomain, mat) result(dir)

        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain
        integer         , intent(in)            :: mat

        !OUTPUT
        integer :: dir

        !LOCAL
        integer :: error, code

        real :: vx, vy, vz

        !!! XXX: this wont work always, read comments from the call site
        vx = Tdomain%sSubDomain(mat)%pml_width(0)
        vy = Tdomain%sSubDomain(mat)%pml_width(1)
        vz = Tdomain%sSubDomain(mat)%pml_width(2)

        !/////////////Defining PML orientation
        dir = 0
        !Face X oriented
        if  (vx /= 0 .and. vy == 0d0 .and. vz == 0d0) then
            dir = 0
        !Face Y oriented
        elseif (vx == 0 .and. vy /= 0d0 .and. vz == 0d0) then
            dir = 1
        !Face Z oriented
        elseif (vx == 0 .and. vy == 0d0 .and. vz /= 0d0) then
            dir = 2
        !Edge in XY
        elseif (vx /= 0 .and. vy /= 0d0 .and. vz == 0d0) then
            dir = 3

        !Edge in YZ
        elseif (vx == 0 .and. vy /= 0d0 .and. vz /= 0d0) then
            dir = 4

        !Edge in ZX
        elseif (vx /= 0 .and. vy == 0d0 .and. vz /= 0d0) then
            dir = 5

        !Vertex in XYZ
        elseif (vx /= 0 .and. vy /= 0d0 .and. vz /= 0d0) then
            dir = 6

        !Undefined PML
        else
            write(*,*) "ERROR in mat ", mat, " (PML) definition (directions), check 'material.input'"
            call MPI_ABORT(Tdomain%communicateur, error, code)
        end if

    end function read_PML_Direction


end module build_prop_files
                     
