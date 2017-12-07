!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module surface_input
    use constants, only : fpp
    implicit none
    
contains
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------

    !!\brief
    !! Carte maîtresse du traitement de surfaces. cette subroutine permet de définir les
    !! les différents paramètres fournit comme jeux de données et notamment l'orientation
    !! des différents problème associée
    !<

    subroutine read_surface_input(Tdomain, config)

        use sdomain
        use ssurf
        use iso_c_binding
        use sem_c_config
        use Alertes

        implicit none
        type(domain), intent(inout)  :: Tdomain
        type(sem_config), intent(in) :: config
        type(sem_surfaces), pointer  :: surf
        character(len=20)            :: char
        character(len=70)            :: sourcename
        real(fpp)                 :: ndir
        integer                      :: i, rg, nsurf, int
        integer, allocatable         :: dummylist(:)
        logical                      :: Boolean(1:20)
        character(len=800)                       :: parametric_var
        character(len=256)           :: FunctionName ='read_surface_input'
        character(len=256)           :: SourceFile = 'create_sem_surface_input'
        character(len=700)           :: ErrorSMS

100     format('         ',A40)

        Tdomain%nsurface = config%nsurface
        allocate (Tdomain%nsurfsource(0:Tdomain%nsurface -1))

        call c_f_pointer(config%surface, surf)
        nsurf          =0
        Tdomain%n_NEBC =0
        Tdomain%n_PWBC =0
        Tdomain%n_FTBC =0
        Tdomain%n_DIRIC=0
        !!
        rg = Tdomain%rank

        if (rg == 0) then
            write(*,*)
            write(*,*) "--> READING SURFACE INPUT"
            write(*,*)
        endif

        do while(associated(surf))

            if ( surf%surface_present /= 0 ) then
                int = nsurf
                write(char,*) int
                if (allocated(dummylist)) deallocate(dummylist)
                if  (surf%surface_whatbc == 1) then
                    Tdomain%nsurfsource(nsurf)%what_bc = 'NE'
                    if (Tdomain%n_NEBC-1.ge.0) allocate(dummylist(1:Tdomain%n_NEBC))
                    if (allocated(Tdomain%list_NEBC)) dummylist = Tdomain%list_NEBC
                    if (allocated(Tdomain%list_NEBC)) deallocate(Tdomain%list_NEBC)
                    Tdomain%n_NEBC = Tdomain%n_NEBC + 1
                    allocate(Tdomain%list_NEBC(1:Tdomain%n_NEBC))
                    if (Tdomain%n_NEBC-1.gt.0) Tdomain%list_NEBC(1:Tdomain%n_NEBC-1) = dummylist
                    Tdomain%list_NEBC(Tdomain%n_NEBC) = nsurf
                    Tdomain%nsurfsource(nsurf)%shape = surf%surface_space
                    Tdomain%nsurfsource(nsurf)%size = surf%surface_size
                    if (rg.eq.0) write(*,100) "Neumann /-> surface"//adjustl(char(1:len_trim(char)))

                elseif (surf%surface_whatbc == 2) then
                    Tdomain%nsurfsource(nsurf)%what_bc = 'PW'
                    if (Tdomain%n_PWBC-1.ge.0) allocate(dummylist(1:Tdomain%n_PWBC))
                    if (allocated(Tdomain%list_PWBC)) dummylist = Tdomain%list_PWBC
                    if (allocated(Tdomain%list_PWBC)) deallocate(Tdomain%list_PWBC)
                    Tdomain%n_PWBC = Tdomain%n_PWBC + 1
                    allocate(Tdomain%list_PWBC(1:Tdomain%n_PWBC))
                    if (Tdomain%n_PWBC-1.gt.0) Tdomain%list_PWBC(1:Tdomain%n_PWBC-1) = dummylist
                    Tdomain%list_PWBC(Tdomain%n_PWBC) = nsurf
                    Tdomain%nsurfsource(nsurf)%wave_type = surf%surface_wave
                    Tdomain%nsurfsource(nsurf)%Speed= surf%surface_Speed
                    Tdomain%nsurfsource(nsurf)%Kdir = surf%surface_dirU
                    if (rg.eq.0) write(*,100) "Plane Wave /-> surface"//adjustl(char(1:len_trim(char)))

                elseif (surf%surface_whatbc == 3) then
                    Tdomain%nsurfsource(nsurf)%what_bc = 'FT'
                    if (Tdomain%n_FTBC-1.ge.0) allocate(dummylist(1:Tdomain%n_FTBC))
                    if (allocated(Tdomain%list_FTBC)) dummylist = Tdomain%list_FTBC
                    if (allocated(Tdomain%list_FTBC)) deallocate(Tdomain%list_FTBC)
                    Tdomain%n_FTBC = Tdomain%n_FTBC + 1
                    if (Tdomain%n_FTBC-1.gt.0) Tdomain%list_FTBC(1:Tdomain%n_FTBC-1) = dummylist
                    allocate(Tdomain%list_FTBC(1:Tdomain%n_FTBC))
                    Tdomain%list_FTBC(Tdomain%n_FTBC) = nsurf
                    if (rg.eq.0) write(*,100)  "Fault /-> surface"//adjustl(char(1:len_trim(char)))

                elseif (surf%surface_whatbc == 4) then
                    Tdomain%nsurfsource(nsurf)%what_bc = 'DR'
                    if (Tdomain%n_DIRIC-1.ge.0) allocate(dummylist(1:Tdomain%n_DIRIC))
                    if (allocated(Tdomain%list_DIRICBC)) dummylist = Tdomain%list_DIRICBC
                    if (allocated(Tdomain%list_DIRICBC)) deallocate(Tdomain%list_DIRICBC)
                    Tdomain%n_DIRIC = Tdomain%n_DIRIC + 1
                    allocate(Tdomain%list_DIRICBC(1:Tdomain%n_DIRIC))
                    if (Tdomain%n_DIRIC-1.gt.0) Tdomain%list_DIRICBC(1:Tdomain%n_DIRIC-1) = dummylist
                    Tdomain%list_DIRICBC(Tdomain%n_DIRIC) = nsurf
                    if (rg.eq.0) write(*,100)  "Dirichlet /-> surface"//adjustl(char(1:len_trim(char)))

                else
                    ErrorSMS = "unknown surface problem (type : ??? ) /-> surface"//adjustl(char(1:len_trim(char)))
                    call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
                endif

                Tdomain%nsurfsource(nsurf)%mat_index = surf%surface_mat

                Boolean = .false.
                do i=lbound(surf%surface_list,1),ubound(surf%surface_list,1)
                    Boolean(i) = surf%surface_list(i) /= 0
                enddo

                ErrorSMS = "unknown associeted surface to "//adjustl(sourcename(1:len_trim(sourcename)))
                if (Count(Boolean) == 0) call ErrorMessage(ErrorSMS,FunctionName,SourceFile)

                allocate(Tdomain%nsurfsource(nsurf)%index(1:Count(Boolean)))
                Tdomain%nsurfsource(nsurf)%index = surf%surface_list(1:Count(Boolean))-1
                if (Tdomain%nsurfsource(nsurf)%what_bc /= 'DR') then
                    if (surf%surface_type==2) then
                        Tdomain%nsurfsource(nsurf)%wtype = 'R'
                        Tdomain%nsurfsource(nsurf)%f0 = surf%surface_f0
                        Tdomain%nsurfsource(nsurf)%Rickertau=surf%Rtau
                        sourcename = "Ricker in time"
                    elseif (surf%surface_type==1) then
                        Tdomain%nsurfsource(nsurf)%wtype = 'G'
                        sourcename = "Gaussian in time"
                    elseif (surf%surface_type==14) then
                        Tdomain%nsurfsource(nsurf)%wtype = 'A'
                    endif

                    ndir = sqrt(surf%surface_L(1)**2+surf%surface_L(2)**2+surf%surface_L(3)**2)
                    if (ndir>0.0) Tdomain%nsurfsource(nsurf)%dir(0:2) = surf%surface_L(1:3)/ndir
                    Tdomain%nsurfsource(nsurf)%scoord(0:2) = surf%surface_C(1:3)
                    Tdomain%nsurfsource(nsurf)%amplitude = surf%amplitude

                    if (Tdomain%nsurfsource(nsurf)%wtype == 'A') then

                        Tdomain%nsurfsource(nsurf)%dim   = surf%surface_dim
                        Tdomain%nsurfsource(nsurf)%varia(1:len_trim(fromcstr(surf%surface_varia))) = trim(fromcstr(surf%surface_varia))
                        Tdomain%nsurfsource(nsurf)%source(1:len_trim(fromcstr(surf%surface_source)))= trim(fromcstr(surf%surface_source))
                        Tdomain%nsurfsource(nsurf)%funcx(1:len_trim(fromcstr(surf%surface_funcx))) = trim(fromcstr(surf%surface_funcx))

                        if (Tdomain%nsurfsource(nsurf)%dim.gt.1) then
                            Tdomain%nsurfsource(nsurf)%funcy(1:len_trim(fromcstr(surf%surface_funcy))) = trim(fromcstr(surf%surface_funcy))
                        endif
                        if (Tdomain%nsurfsource(nsurf)%dim.gt.2) then
                            Tdomain%nsurfsource(nsurf)%funcz(1:len_trim(fromcstr(surf%surface_funcz))) = trim(fromcstr(surf%surface_funcz))
                        endif
                        if (Tdomain%nsurfsource(nsurf)%source == 'M') then
                            Tdomain%nsurfsource(nsurf)%funcxy(1:len_trim(fromcstr(surf%surface_funcxy)))= trim(fromcstr(surf%surface_funcxy))
                            if (Tdomain%nsurfsource(nsurf)%dim == 3) then
                                Tdomain%nsurfsource(nsurf)%funcxz(1:len_trim(fromcstr(surf%surface_funcxz)))= trim(fromcstr(surf%surface_funcxz))
                                Tdomain%nsurfsource(nsurf)%funcyz(1:len_trim(fromcstr(surf%surface_funcyz)))= trim(fromcstr(surf%surface_funcyz))
                            endif
                        endif
                        !!
                        !!
                        Tdomain%nsurfsource(nsurf)%paramvar = surf%surface_paramvar
                        if (Tdomain%nsurfsource(nsurf)%paramvar==1) then
                            Tdomain%nsurfsource(nsurf)%nparamvar = surf%surface_nparamvar
                            allocate(Tdomain%nsurfsource(nsurf)%paravalue(1:surf%surface_nparamvar))
                            do i=1,surf%surface_nparamvar
                                Tdomain%nsurfsource(nsurf)%paravalue(i) = surf%surface_Paravalue(i)
                            enddo
                            parametric_var(1:len_trim(fromcstr(surf%surface_paramname))) = trim(fromcstr(surf%surface_paramname))
                            call Split(parametric_var, surf%surface_nparamvar, Tdomain%nsurfsource(nsurf)%paramname," ")
                        endif
                        !!
                        !!
                        if (rg.eq.0) then
                            write(*,*)
                            write(*,*) "Surface BC Source analytical defined. "
                            write(*,*) "-------------------------------------"
                            write(*,'(A12)') "Va = "// adjustr(Tdomain%nsurfsource(nsurf)%varia)
                            write(*,'(A150)') "F1 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcx)
                            if (Tdomain%nsurfsource(nsurf)%dim.gt.1) then
                                write(*,'(A150)') "F2 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcy)
                            endif
                            if (Tdomain%nsurfsource(nsurf)%dim.gt.2) then
                                write(*,'(A150)') "F3 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcz)
                            endif
                            if (Tdomain%nsurfsource(nsurf)%source.eq."M") then
                                write(*,'(A150)') "F4 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcxy)
                                if (surf%surface_dim.gt.2) then
                                    write(*,'(A150)') "F5 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcxz)
                                    write(*,'(A150)') "F6 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcyz)
                                endif
                            endif
                            if (Tdomain%nsurfsource(nsurf)%paramvar==1) then
                                write(*,*)
                                write(*,'(2A12)') "Param"," ParamVal"
                                do i=1,Tdomain%nsurfsource(nsurf)%nparamvar
                                    write(*,'(         A2     ES13.6)') Tdomain%nsurfsource(nsurf)%paramname(i), &
                                        Tdomain%nsurfsource(nsurf)%paravalue(i)
                                enddo
                            endif
                        endif
                        write(*,*)

                    else

                        if (rg==0) then
                            write(*,*) "Surface BC Source  : "//adjustl(sourcename(1:len_trim(sourcename)))
                            write(*,'(A22,3ES13.3)') " Direction          : ", Tdomain%nsurfsource(nsurf)%dir

                            sourcename = fromcstr(surf%surface_name)
                            write(*,'(A22,A12,1ES13.3,A1)') " Load shape         : ", sourcename(1:len_trim(sourcename))// &
                                adjustl(" (R="),Tdomain%nsurfsource(nsurf)%size,")"
                            if (Tdomain%nsurfsource(nsurf)%f0/=0) write(*,'(A22,ES13.3)') "Ricker Freqence    : ", &
                                Tdomain%nsurfsource(nsurf)%f0
                            if (Tdomain%nsurfsource(nsurf)%Rickertau/=0) write(*,'(A22,ES13.3)') " Ricker tau         : ", &
                                Tdomain%nsurfsource(nsurf)%Rickertau
                        endif
                    endif

                    if (rg==0) then
                        write(*,'(A21,3ES13.6)') " Source Ref coord   : ", Tdomain%nsurfsource(nsurf)%scoord
                        write(*,'(A21,30I4)') " Associated surface : ",(Tdomain%nsurfsource(nsurf)%index(i), i=1,Count(Boolean))
                        write(*,'(A21,ES13.6)') " Load amplitude     : ", Tdomain%nsurfsource(nsurf)%amplitude
                    endif

                else
                    write(*,'(A21,30I4)') "Associated surface : ",(Tdomain%nsurfsource(nsurf)%index(i), i=1,Count(Boolean))
                    write(*,*)
                endif
                nsurf = nsurf + 1
                Tdomain%logicD%surfBC = .true.
                write(*,*)
            endif
            call c_f_pointer(surf%next, surf)

        enddo

        if (Tdomain%nsurface /= nsurf) Tdomain%nsurface = nsurf
        if (allocated(dummylist)) deallocate(dummylist)
        if (nsurf==0) Tdomain%logicD%surfBC = .false.
        write(*,*)

    end subroutine read_surface_input
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine Split(String, nstring, SubStrings,Delimiter)

        implicit none

        character(len=*), intent(in)                         :: String
        integer, intent(in)                                  :: nstring
        character, intent(in)                                :: Delimiter
        character(len=*), dimension(nstring), intent(inout)  :: SubStrings
        integer                                              :: i, j, k, l, n, ns
        logical                                              :: EndOfLine

        do i=1,nstring
            SubStrings(i) = ' '
        end do
        n = 0
        k = 1
        l = len_trim(String)
        EndOfLine = l-k < 0
        do while (.not.EndOfLine)
            j = index(String(k:l),Delimiter)
            if (j == 0) then
                j=l+1
            else
                j=j+k-1
            end if
            n=n+1
            if ((j.ne.k).and.(len_trim(String(k:j-1)).ne.0).and.(n.le.nstring)) then
                SubStrings(n) = String(k:j-1)
                ns =n;
            else
                n  =ns;
            endif
            k=j+1
            EndOfLine=l-k < 0
        end do

    end subroutine Split
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------

    !!\brief
    !! Subroutine permettant d'associée au différentes surfaces les propriétés élastiques en fonction
    !! des materiaux sur lesquel elles sont affectées
    !<

    subroutine define_surface_properties(Tdomain)

        use sdomain
        use Alertes

        implicit none
        type(domain), intent(inout) :: Tdomain
        integer                     :: ns, index, surf , s, n
        character(len=20)           :: char
        character(len=256)          :: FunctionName ='define_surface_properties'
        character(len=256)          :: SourceFile = 'create_sem_surface_input'
        character(len=700)          :: ErrorSMS

2004    format( A20,I4,A4,A)
2007    format('        |',A12, I6)
2009    format('        |',A12, ES13.6)

        ! elastic properties
        do ns = 0, Tdomain%nsurface-1
            index = Tdomain%nsurfsource(ns)%mat_index
            do s = 1,size(Tdomain%nsurfsource(ns)%index)
                surf = Tdomain%nsurfsource(ns)%index(s)
                write(char,*) surf
                do n=0,size(Tdomain%sSurfaces)-1
                    if (Tdomain%sSurfaces(n)%name == "surface"//adjustl(char(:len_trim(char)))) then
                        Tdomain%sSurfaces(n)%Elastic%lambda  = Tdomain%sSubDomain(index)%DLambda
                        Tdomain%sSurfaces(n)%Elastic%Mu      = Tdomain%sSubDomain(index)%DMu
                        Tdomain%sSurfaces(n)%Elastic%Sspeed  = Tdomain%sSubDomain(index)%Sspeed
                        Tdomain%sSurfaces(n)%Elastic%Pspeed  = Tdomain%sSubDomain(index)%Pspeed
                        Tdomain%sSurfaces(n)%Elastic%density = Tdomain%sSubDomain(index)%dDensity
                        if (Tdomain%nsurfsource(ns)%what_bc=='PW') then
                            call PlaneWaveSpeed(Tdomain%nsurfsource(ns),Tdomain%sSurfaces(n))
                        endif
                        !! print parameters values in the screen
                        ErrorSMS = "unknown associeted domain to "//adjustl(Tdomain%sSurfaces(n)%name)
                        if (Tdomain%nsurfsource(ns)%mat_index<0) call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
                        if (Tdomain%rank == 0) then
                            write(*,*)
                            write(*,2004) "SURFACE -> ", surf
                            write(*,2007) "Domaim  : " , Tdomain%nsurfsource(ns)%mat_index
                            write(*,2009) "Lambda  : " , Tdomain%sSurfaces(n)%Elastic%lambda
                            write(*,2009) "Mu      : " , Tdomain%sSurfaces(n)%Elastic%Mu
                            write(*,2009) "Speed Vp: ",  Tdomain%sSurfaces(n)%Elastic%Pspeed
                            write(*,2009) "Speed Vs: ",  Tdomain%sSurfaces(n)%Elastic%Sspeed
                            if (Tdomain%nsurfsource(ns)%what_bc=='PW') then
                                write(*,2009) "Speed PW: ", Tdomain%sSurfaces(n)%Elastic%PWspeed
                                write(*,2009) "PW dirK : ", Tdomain%nsurfsource(ns)%dir
                                write(*,2009) "PW dirD : ", Tdomain%nsurfsource(ns)%Kdir
                            endif
                            write(*,2009) "Density : ",  Tdomain%sSurfaces(n)%Elastic%density
                        endif
                    endif
                enddo
            enddo
        end do
        write(*,*)

    end subroutine define_surface_properties
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------

    !!\brief
    !! Subroutine qui permet de définir les principale directions des ondes plane classiques
    !! et aussi les différentes vitesses de propagations qui leur sont associées
    !<

    subroutine PlaneWaveSpeed(surfsrc,surf)

        use ssurf

        implicit none
        type(SurfaceParam),        intent(inout) :: surfsrc
        type(SurfaceT),            intent(inout) :: surf
        real(fpp), dimension(0:2)             :: dirP, dirD

        dirP = surfsrc%dir

        if (surfsrc%wave_type==1) then
            dirD = dirP
            surf%Elastic%PWspeed = surf%Elastic%Pspeed
            write(*,*)
            write(*,*) " LONGITUDINAL INCIDENT PLANE WAVE "
            write(*,*)
        elseif (surfsrc%wave_type==2) then
            dirD = DirP
            if ((dirP(1)**2 + dirP(2)**2) /=0.d0) dirD(0) = -dirP(0) /(dirP(1)**2 + dirP(2)**2 )
            if ((dirP(0)**2 + dirP(2)**2) /=0.d0) dirD(0) = -dirP(1) /(dirP(0)**2 + dirP(2)**2 )
            if ((dirP(0)**2 + dirP(1)**2) /=0.d0) dirD(0) = -dirP(0) /(dirP(0)**2 + dirP(1)**2 )
            surf%Elastic%PWspeed = surf%Elastic%Sspeed
            write(*,*)
            write(*,*) " TRANSVERSAL INCIDENT PLANE WAVE "
            write(*,*)
        elseif (surfsrc%wave_type==3) then
            dirD(0)=0.d0; dirD(1)=0.d0; dirD(2) = 1.d0
            surf%Elastic%PWspeed = surf%Elastic%Sspeed
            write(*,*)
            write(*,*) " SH-WAVE INCIDENT PLANE WAVE "
            write(*,*)
        elseif (surfsrc%wave_type==4) then
            dirD(0)=-dirP(1); dirD(1)=dirP(0); dirD(2)=0.d0
            surf%Elastic%PWspeed = surf%Elastic%Sspeed
            write(*,*)
            write(*,*) " SV-WAVE INCIDENT PLANE WAVE "
            write(*,*)
        else
            surf%Elastic%PWspeed = surf%Elastic%Sspeed
            dirD = surfsrc%Kdir
            write(*,*)
            write(*,*) " SV-WAVE INCIDENT PLANE WAVE "
            write(*,*)
        endif

        surfsrc%Kdir =  dirD

    end subroutine PlaneWaveSpeed
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------

    !!\brief
    !! Subroutine permettant de vérifier si les tags de surfaces fournit par l'utilisateur sont
    !! sont les même que ceux lui dans les mesh4spec. Cette routine nécessite d'être paralléliser
    !! car toutes les surfaces ne figures pas surtous les proc
    !<

    subroutine  surface_in_list(Tdomain)

        use sdomain
        use mpi
        use Alertes

        implicit none
        type(domain), intent(in)                   :: Tdomain
        integer                                    :: n, ns, s, intt, code, rg
        integer                                    :: nb_proc
        logical                                    :: find
        logical,dimension(:), allocatable          :: Tab_logic
        integer ,dimension(:), allocatable         :: surf_tags
        integer, dimension(:), allocatable         :: Received, Depla
        integer, dimension(:), allocatable         :: list_index
        character(len=100)          :: surflist_index
        character(len=256)          :: FunctionName ='surface_in_list'
        character(len=256)          :: SourceFile = 'create_sem_surface_input'
        character(len=700)          :: ErrorSMS
        character(len=20)           :: charr
        integer                     :: surfnum

        !! Blocks until all processes have reached this step
        call MPI_Barrier(Tdomain%communicateur, code)
        call MPI_Comm_size(MPI_COMM_WORLD, nb_proc,code)
        call MPI_Comm_Rank(MPI_COMM_WORLD, rg, code)

        allocate(Received(1:nb_proc), Depla(1:nb_proc))
        Received = 0
        ! On recherchele nombre de surface lu sur chaque proc

        call MPI_ALLGATHER(size(Tdomain%sSurfaces),1,MPI_INTEGER,Received(:),1, MPI_INTEGER,MPI_COMM_WORLD,code)

        if (rg==0) then
            Depla(:)=0
            do n=2,nb_proc
                Depla(n)=sum(Received(1:n-1))
            enddo
        endif
        allocate(Tab_logic(nb_proc), surf_tags(sum(Received)), list_index(size(Tdomain%sSurfaces)))
        list_index =0; Tab_logic = .true.; surf_tags =0

        do ns = 0, Tdomain%nsurface-1
            do s = 1,size(Tdomain%nsurfsource(ns)%index)
                find = .false.
                intt = Tdomain%nsurfsource(ns)%index(s)
                write(charr,'(I16)') intt
                do n=0,size(Tdomain%sSurfaces)-1
                    read(Tdomain%sSurfaces(n)%name(8:len_trim(Tdomain%sSurfaces(n)%name)),*) surfnum
                    list_index(n+1) = surfnum
                    if (Tdomain%sSurfaces(n)%name == "surface"//adjustl(charr(1:len_trim(charr)))) then
                        find = .true.
                        exit
                    endif
                enddo
                ! Récuperation et diffussion de la valeur de find sur tous les proc
                call MPI_ALLGATHER(find,1,MPI_LOGICAL,Tab_logic,1, MPI_LOGICAL,MPI_COMM_WORLD,code)
                ! liste complète des tags retransmise sur le proc 0, surf_tags va contenir toutes les listes
                ! récupérées sur tous les proc. cette liste peut contenir des doublons
                call MPI_GATHERV(list_index, size(list_index), MPI_INTEGER, surf_tags, Received, Depla, MPI_INTEGER, 0, MPI_COMM_WORLD, code)

                if ((rg == 0).and.(count(Tab_logic)==0)) then
                    call Unique_string(surf_tags,surflist_index)
                    intt = intt+1; write(charr,'(I16)') intt
                    ErrorSMS = adjustr("surface (index ="//adjustl(charr))//adjustl(") is not found in tag list: ")// &
                        adjustl(surflist_index(:len_trim(surflist_index)))
                    call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
                endif
            enddo
        enddo
        deallocate(Received, Depla, Tab_logic, surf_tags, list_index)

    end subroutine surface_in_list
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------

    !!\brief
    !! Reconstruire la liste de tags de surfaces retrouvés dans mesh4spec tout en évitant les
    !! les répétitions. string est un tableau de chaine de caratères et de taille nb_proc.
    !! Chaque composante de ce tableau est une liste des tags de surfaces retrouvés sur chaque
    !! proc. Cette fonction a la même fonctionnalité que celle de C++ (std::unique)
    !<

    subroutine  Unique_string(string,Unique)

        implicit none

        integer, dimension(:),             intent(in) :: string
        character(len=*) , intent(out)                :: Unique
        character(len=10), dimension(:), allocatable :: res
        character(len=150), dimension(:), allocatable :: Character_list
        character(len=20)                             :: char
        integer                                       :: k, i, j

        write(char,*) string(1)
        unique = char
        do i=2,size(string)
            write(char,*) string(i)
            Unique=adjustr(Unique(1:len_trim(Unique))//";")//adjustl(char(1:len_trim(char)))
        enddo
        j=0
        do i=1,len_trim(Unique)
            if (Unique(i:i)==";") j=j+1
        enddo
        allocate(Character_list(j))
        call Split(Unique,j,Character_list,";")
        k=1;
        allocate(res(1:i))
        res(1)=Character_list(1)
        Unique=adjustl(res(1))
        block : &
            do i=2,size(Character_list)
          do j=1,k
              if ((trim(res(j))==trim(Character_list(i))).or.(len_trim(Character_list(i))==0)) then
                  cycle block
              endif
          enddo
          Unique=adjustr(Unique(1:len_trim(Unique))//";")//adjustl(Character_list(i))
          k=k+1
          res(k)=Character_list(i)
      enddo block
      deallocate(res,Character_list)

  end subroutine Unique_string

end module surface_input

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
