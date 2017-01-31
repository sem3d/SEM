!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file save_deformation.F90
!!\brief Contient la subroutine save_deformation.
!!
!! Calcul de la deformation
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<
#ifdef MKA3D
subroutine save_deformation (Tdomain,it,i_snap,sortie_capteur_deformation,nom_dir_sorties)
#else
subroutine save_deformation (Tdomain,it,i_snap,sortie_capteur_deformation)
#endif
    use sdomain
    use semdatafiles
    use mpi

    implicit none
    type (domain), intent (INOUT):: TDomain
    integer, intent (IN) :: it,i_snap
    logical :: sortie_capteur_deformation

    ! local variables
    integer :: ngllx,ngllz,ipoint, i,j,n, mat,nf
    real(fpp) :: epsilon_plus,epsilon_minus,zzero
    real(fpp), dimension(:,:) ,allocatable  :: Grad_u
    real(fpp), dimension (:), allocatable  :: Wheight
    real(fpp), dimension (:,:), allocatable :: Disp
    real(fpp), dimension (:,:,:), allocatable :: Loc_grad_U , Grad_x, Displ
    real(fpp), dimension (:,:), pointer :: hTprime, hprime

    ! Variable from MPI
    integer :: i_proc, n_face_pointed, ne , wf, i_send, i_stock, tag_send, tag_receive, ierr
    integer , dimension  (MPI_STATUS_SIZE) :: status
    real(fpp), dimension (0:4, 0:10000) :: Aus
    real(fpp), dimension (:,:), allocatable :: Send_data, Receive_data

    !modif CSSI: correction bug
    integer,dimension(:),allocatable :: contributeur

    character(len=MAX_FILE_SIZE) :: fnamef

#ifdef MKA3D
    character(len=24), intent(IN) :: nom_dir_sorties
#endif




    allocate (Grad_u(0:3,0:Tdomain%n_glob_points-1))
    allocate (Wheight(0:Tdomain%n_glob_points-1))
    zzero = 0. !!Ajout Gsa sinon erreur a l'execution parfois
    Grad_u = 0
    Wheight = 0
    do n = 0,Tdomain%n_elem-1
        if (.not. Tdomain%specel(n)%PML) then
            ngllx = Tdomain%specel(n)%ngllx;  ngllz = Tdomain%specel(n)%ngllz
            allocate (Loc_Grad_U(0:ngllx-1,0:ngllz-1,0:3))
            allocate (Grad_x(0:ngllx-1,0:ngllz-1,0:3))
            allocate  (Disp (0:ngllx-1,0:ngllz-1))
            allocate (Displ(0:ngllx-1,0:ngllz-1,0:1))

            nf = Tdomain%specel(n)%near_face(0)

            if ( Tdomain%sFace(nf)%near_element(0) == n) then
                Displ(1:ngllx-2,0,0) =  Tdomain%sFace(nf)%Displ(:,0)
                Displ(1:ngllx-2,0,1) =  Tdomain%sFace(nf)%Displ(:,1)
            else
                if ( Tdomain%sFace(nf)%coherency) then
                    Displ(1:ngllx-2,0,0) =  Tdomain%sFace(nf)%Displ(:,0)
                    Displ(1:ngllx-2,0,1) =  Tdomain%sFace(nf)%Displ(:,1)
                else
                    do i = 1, ngllx-2
                        Displ(i,0,0) =  Tdomain%sFace(nf)%Displ(ngllx-1-i,0)
                        Displ(i,0,1) =  Tdomain%sFace(nf)%Displ(ngllx-1-i,1)
                    enddo
                endif
            endif


            nf = Tdomain%specel(n)%near_face(1)

            if ( Tdomain%sFace(nf)%near_element(0) == n) then
                Displ(ngllx-1,1:ngllz-2,0) =  Tdomain%sFace(nf)%Displ(:,0)
                Displ(ngllx-1,1:ngllz-2,1) =  Tdomain%sFace(nf)%Displ(:,1)
            else
                if ( Tdomain%sFace(nf)%coherency) then
                    Displ(ngllx-1,1:ngllz-2,0)=  Tdomain%sFace(nf)%Displ(:,0)
                    Displ(ngllx-1,1:ngllz-2,1)=  Tdomain%sFace(nf)%Displ(:,1)
                else
                    do i = 1, ngllz-2
                        Displ(ngllx-1,i,0) =  Tdomain%sFace(nf)%Displ(ngllz-1-i,0)
                        Displ(ngllx-1,i,1) =  Tdomain%sFace(nf)%Displ(ngllz-1-i,1)
                    enddo
                endif
            endif


            nf = Tdomain%specel(n)%near_face(2)

            if ( Tdomain%sFace(nf)%near_element(0) == n) then
                Displ(1:ngllx-2,ngllz-1,0) =  Tdomain%sFace(nf)%Displ(:,0)
                Displ(1:ngllx-2,ngllz-1,1) =  Tdomain%sFace(nf)%Displ(:,1)
            else
                if ( Tdomain%sFace(nf)%coherency) then
                    Displ(1:ngllx-2,ngllz-1,0) =  Tdomain%sFace(nf)%Displ(:,0)
                    Displ(1:ngllx-2,ngllz-1,1) =  Tdomain%sFace(nf)%Displ(:,1)
                else
                    do i = 1, ngllx-2
                        Displ(i,ngllz-1,0) =  Tdomain%sFace(nf)%Displ(ngllx-1-i,0)
                        Displ(i,ngllz-1,1) =  Tdomain%sFace(nf)%Displ(ngllx-1-i,1)
                    enddo
                endif
            endif


            nf = Tdomain%specel(n)%near_face(3)

            if ( Tdomain%sFace(nf)%near_element(0) == n) then
                Displ(0,1:ngllz-2,0) =  Tdomain%sFace(nf)%Displ(:,0)
                Displ(0,1:ngllz-2,1) =  Tdomain%sFace(nf)%Displ(:,1)
            else
                if ( Tdomain%sFace(nf)%coherency) then
                    Displ(0,1:ngllz-2,0)=  Tdomain%sFace(nf)%Displ(:,0)
                    Displ(0,1:ngllz-2,1)=  Tdomain%sFace(nf)%Displ(:,1)
                else
                    do i = 1, ngllz-2
                        Displ(0,i,0) =  Tdomain%sFace(nf)%Displ(ngllz-1-i,0)
                        Displ(0,i,1) =  Tdomain%sFace(nf)%Displ(ngllz-1-i,1)
                    enddo
                endif
            endif





            nf = Tdomain%specel(n)%Near_Vertex(0)
            Displ(0,0,0:1) = Tdomain%sVertex(nf)%Displ(0:1)

            nf = Tdomain%specel(n)%Near_Vertex(1)
            Displ(ngllx-1,0,0:1) = Tdomain%sVertex(nf)%Displ(0:1)

            nf = Tdomain%specel(n)%Near_Vertex(2)
            Displ(ngllx-1,ngllz-1,0:1) = Tdomain%sVertex(nf)%Displ(0:1)

            nf = Tdomain%specel(n)%Near_Vertex(3)
            Displ(0,ngllz-1,0:1) = Tdomain%sVertex(nf)%Displ(0:1)

            Displ(1:ngllx-2,1:ngllz-2,0:1) = Tdomain%specel(n)%Displ(:,:,0:1)

            mat = Tdomain%specel(n)%mat_index
            hTprime => Tdomain%sSubdomain(mat)%HTprimex
            hprime => Tdomain%sSubdomain(mat)%Hprimez

            Disp=Displ(:,:,0)
            Loc_Grad_U(:,:,0) =MATMUL (hTprime, Disp)   ! dux_dxi
            Loc_Grad_U (:,:,1) = MATMUL (Disp, hprime) ! dux_deta
            Disp = Displ(:,:,1)
            Loc_Grad_U(:,:,2) =MATMUL (hTprime, Disp)   ! duz_dxi
            Loc_Grad_U (:,:,3) = MATMUL (Disp, hprime) ! duz_deta

            Grad_x (:,:,0) = Loc_Grad_U(:,:,0) * Tdomain%specel(n)%Invgrad(:,:,0,0) + Loc_Grad_U (:,:,1) *  &
                Tdomain%specel(n)%Invgrad(:,:,0,1)   ! dux_dx
            Grad_x (:,:,1) = Loc_Grad_U(:,:,0) * Tdomain%specel(n)%Invgrad(:,:,1,0) + Loc_Grad_U (:,:,1) *  &
                Tdomain%specel(n)%Invgrad(:,:,1,1)   ! dux_dz
            Grad_x (:,:,2) = Loc_Grad_U(:,:,2) * Tdomain%specel(n)%Invgrad(:,:,0,0) + Loc_Grad_U (:,:,3) *  &
                Tdomain%specel(n)%Invgrad(:,:,0,1)  ! duz_dx
            Grad_x (:,:,3) = Loc_Grad_U(:,:,2) * Tdomain%specel(n)%Invgrad(:,:,1,0) + Loc_Grad_U (:,:,3) *  &
                Tdomain%specel(n)%Invgrad(:,:,1,1)! duz_dz


            do j = 0,ngllz-1
                do i = 0,ngllx-1
                    ipoint = Tdomain%specel(n)%Iglobnum(i,j)

                    Grad_u (0,ipoint) = Grad_x(i,j,0) * Tdomain%specel(n)%Jacob(i,j) +   Grad_u (0,ipoint)
                    Grad_u (1,ipoint) = Grad_x(i,j,1) * Tdomain%specel(n)%Jacob(i,j) +   Grad_u (1,ipoint)
                    Grad_u (2,ipoint) = Grad_x(i,j,2) * Tdomain%specel(n)%Jacob(i,j) +   Grad_u (2,ipoint)
                    Grad_u (3,ipoint) = Grad_x(i,j,3) * Tdomain%specel(n)%Jacob(i,j) +   Grad_u (3,ipoint)
                    Wheight (ipoint) = Wheight (ipoint) + Tdomain%specel(n)%Jacob(i,j)

                enddo
            enddo
            deallocate (Loc_Grad_U, Grad_x, Disp, Displ)
        endif
    enddo

    !modif CSSI: correction bug pour les pt de gauss sur la frontiere du domaine commune avec un autre proc:
    !Grad_u est calcule en chaque pt de gauss sur le pro courant. on boucle sur les faces wall (communes a 2 proc) pour remplir le tableau d envoi aus qui va contenir
    !la contribution Gradu et weight des pt de gauss qui sont sur les faces wall. Mais attention,on ne doit collecter qu'un seule contrib (gradu,w) par pt de gauss !
    !Or qd on boucle sur les faces, on passe plusieurs fois sur le pdgauss, il faut donc tester si un pdg a deja contrinuer avant de le stocker => tableau contributeur

    allocate(contributeur(0:Tdomain%n_glob_points-1))


    ! Communicating values among processors
    do i_proc = 0, Tdomain%n_communications - 1
        i_send = Tdomain%Communication_list (i_proc)
        i_stock = 0
        contributeur(:) = 0
        do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
            n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
            ne = Tdomain%sFace(n_face_pointed)%Near_Element(0)
            wf =Tdomain%sFace(n_face_pointed)%Which_face(0)
            ngllx = Tdomain%specel(ne)%ngllx
            ngllz = Tdomain%specel(ne)%ngllz
            select case (wf)
            case (0)
                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    do i =0, ngllx-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(i,0)
                        Aus (0:3,i_stock) = Grad_u(0:3,ipoint)
                        !Aus (4,i_stock) = Tdomain%specel(ne)%Jacob(i,0)
                        Aus (4,i_stock) = Wheight (ipoint)
                        contributeur(ipoint)=contributeur(ipoint)+1
                        if ( contributeur(ipoint).gt.1 ) Aus (0:4,i_stock) = 0.
                        i_stock = i_stock+1
                    enddo
                else
                    do i = ngllx-1,0,-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(i,0)
                        Aus (0:3,i_stock) = Grad_u(0:3,ipoint)
                        !Aus (4,i_stock) = Tdomain%specel(ne)%Jacob(i,0)
                        Aus (4,i_stock) = Wheight (ipoint)
                        contributeur(ipoint)=contributeur(ipoint)+1
                        if ( contributeur(ipoint).gt.1 ) Aus (0:4,i_stock) = 0.
                        i_stock = i_stock+1
                    enddo
                endif
            case (2)
                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    do i =0, ngllx-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(i,ngllz-1)
                        Aus (0:3,i_stock) = Grad_u(0:3,ipoint)
                        !Aus (4,i_stock) = Tdomain%specel(ne)%Jacob(i,ngllz-1)
                        Aus (4,i_stock) = Wheight (ipoint)
                        contributeur(ipoint)=contributeur(ipoint)+1
                        if ( contributeur(ipoint).gt.1 ) Aus (0:4,i_stock) = 0.
                        i_stock = i_stock+1
                    enddo
                else
                    do i = ngllx-1,0,-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(i,ngllz-1)
                        Aus (0:3,i_stock) = Grad_u(0:3,ipoint)
                        !Aus (4,i_stock) = Tdomain%specel(ne)%Jacob(i,ngllz-1)
                        Aus (4,i_stock) = Wheight (ipoint)
                        contributeur(ipoint)=contributeur(ipoint)+1
                        if ( contributeur(ipoint).gt.1 ) Aus (0:4,i_stock) = 0.
                        i_stock = i_stock+1
                    enddo
                endif
            case (1)
                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    do i =0, ngllz-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(ngllx-1,i)
                        Aus (0:3,i_stock) = Grad_u(0:3,ipoint)
                        !Aus (4,i_stock) = Tdomain%specel(ne)%Jacob(ngllx-1,i)
                        Aus (4,i_stock) = Wheight (ipoint)
                        contributeur(ipoint)=contributeur(ipoint)+1
                        if ( contributeur(ipoint).gt.1 ) Aus (0:4,i_stock) = 0.
                        i_stock = i_stock+1
                    enddo
                else
                    do i = ngllz-1,0,-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(ngllx-1,i)
                        Aus (0:3,i_stock) = Grad_u(0:3,ipoint)
                        !Aus (4,i_stock) = Tdomain%specel(ne)%Jacob(ngllx-1,i)
                        Aus (4,i_stock) = Wheight (ipoint)
                        contributeur(ipoint)=contributeur(ipoint)+1
                        if ( contributeur(ipoint).gt.1 ) Aus (0:4,i_stock) = 0.
                        i_stock = i_stock+1
                    enddo
                endif
            case (3)
                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    do i =0, ngllz-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(0,i)
                        Aus (0:3,i_stock) = Grad_u(0:3,ipoint)
                        !Aus (4,i_stock) = Tdomain%specel(ne)%Jacob(0,i)
                        Aus (4,i_stock) = Wheight (ipoint)
                        contributeur(ipoint)=contributeur(ipoint)+1
                        if ( contributeur(ipoint).gt.1 ) Aus (0:4,i_stock) = 0.

                        i_stock = i_stock+1

                    enddo
                else
                    do i = ngllz-1,0,-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(0,i)
                        Aus (0:3,i_stock) = Grad_u(0:3,ipoint)
                        !Aus (4,i_stock) = Tdomain%specel(ne)%Jacob(0,i)
                        Aus (4,i_stock) = Wheight (ipoint)
                        contributeur(ipoint)=contributeur(ipoint)+1
                        if ( contributeur(ipoint).gt.1 ) Aus (0:4,i_stock) = 0.
                        i_stock = i_stock+1
                    enddo
                endif
            end select
        enddo


        allocate (Send_data(0:2,0:i_stock-1))
        allocate (Receive_data(0:2,0:i_stock-1))
        Send_data (0:2,0:i_stock-1) = Aus (0:2,0:i_stock-1)
        tag_send = i_send * Tdomain%MPI_var%n_proc +Tdomain%MPI_var%my_rank + 100
        tag_receive = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 100

        call MPI_SENDRECV(Send_data, 3*i_stock,MPI_DOUBLE_PRECISION,i_send, tag_send, Receive_data, &
            3*i_stock, MPI_DOUBLE_PRECISION,i_send,tag_receive,Tdomain%communicateur, status, ierr )

        Aus (0:2,0:i_stock-1) = Receive_data(:,:)
        Send_data = 0
        Receive_data = 0
        Send_data (0:1,0:i_stock-1) = Aus (3:4,0:i_stock-1)

        tag_send = i_send * Tdomain%MPI_var%n_proc +Tdomain%MPI_var%my_rank +200
        tag_receive = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 200

        call MPI_SENDRECV(Send_data, 3*i_stock,MPI_DOUBLE_PRECISION,i_send, tag_send, Receive_data, &
            3*i_stock, MPI_DOUBLE_PRECISION,i_send,tag_receive,Tdomain%communicateur, status, ierr )

        Aus (3:4,0:i_stock-1) = Receive_data(0:1,0:i_stock-1)
        deallocate (Receive_data)
        allocate (Receive_data(0:4,0:i_stock-1))
        Receive_data(0:4,0:i_stock-1) = Aus (0:4,0:i_stock-1)

        i_stock = 0
        do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
            n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
            ne = Tdomain%sFace(n_face_pointed)%Near_Element(0)
            wf =Tdomain%sFace(n_face_pointed)%Which_face(0)
            ngllx = Tdomain%specel(ne)%ngllx
            ngllz = Tdomain%specel(ne)%ngllz
            select case (wf)
            case (0)
                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    do i =0, ngllx-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(i,0)

                        Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock)
                        ! Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock) * Receive_data (4,i_stock)
                        Wheight (ipoint) = Wheight (ipoint) + Receive_data (4,i_stock)
                        i_stock = i_stock+1
                    enddo
                else
                    do i = ngllx-1,0,-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(i,0)
                        Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock)
                        ! Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock) * Receive_data (4,i_stock)
                        Wheight (ipoint) = Wheight (ipoint) + Receive_data (4,i_stock)
                        i_stock = i_stock+1
                    enddo
                endif
            case (2)
                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    do i =0, ngllx-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(i,ngllz-1)

                        Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock)
                        ! Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock) * Receive_data (4,i_stock)
                        Wheight (ipoint) = Wheight (ipoint) + Receive_data (4,i_stock)

                        i_stock = i_stock+1
                    enddo
                else
                    do i = ngllx-1,0,-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(i,ngllz-1)
                        Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock)
                        !  Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock) * Receive_data (4,i_stock)
                        Wheight (ipoint) = Wheight (ipoint) + Receive_data (4,i_stock)
                        i_stock = i_stock+1
                    enddo
                endif
            case (1)
                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    do i =0, ngllz-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(ngllx-1,i)
                        Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock)
                        !Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock) * Receive_data (4,i_stock)
                        Wheight (ipoint) = Wheight (ipoint) + Receive_data (4,i_stock)

                        i_stock = i_stock+1
                    enddo
                else
                    do i = ngllz-1,0,-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(ngllx-1,i)
                        Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock)
                        !Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock) * Receive_data (4,i_stock)
                        Wheight (ipoint) = Wheight (ipoint) + Receive_data (4,i_stock)
                        i_stock = i_stock+1
                    enddo
                endif
            case (3)
                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    do i =0, ngllz-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(0,i)
                        Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock)
                        !Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock) * Receive_data (4,i_stock)
                        Wheight (ipoint) = Wheight (ipoint) + Receive_data (4,i_stock)

                        i_stock = i_stock+1
                    enddo
                else
                    do i = ngllz-1,0,-1
                        ipoint = Tdomain%specel(ne)%Iglobnum(0,i)
                        Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock)
                        !Grad_u(0:3,ipoint) = Grad_u(0:3,ipoint) + Receive_data(0:3,i_stock) * Receive_data (4,i_stock)
                        Wheight (ipoint) = Wheight (ipoint) + Receive_data (4,i_stock)
                        i_stock = i_stock+1
                    enddo
                endif
            end select
        enddo
        deallocate (Send_data)
        deallocate (Receive_data)
    enddo

    deallocate(contributeur)





    do n = 0, Tdomain%n_glob_points-1
        if (abs(Grad_u(0,n)) < 1e-40) Grad_u(0,n) = 0.
        if (abs(Grad_u(1,n)) < 1e-40) Grad_u(1,n) = 0.
        if (abs(Grad_u(2,n)) < 1e-40) Grad_u(2,n) = 0.
        if (abs(Grad_u(3,n)) < 1e-40) Grad_u(3,n) = 0.

        if (Wheight (n) > 0 )  Grad_u(:,n) = Grad_u(:,n) / Wheight(n)
    enddo

    if (sortie_capteur_deformation) then

        ! AC initialisation des deformations pour les capteurs
        Tdomain%GrandeurDeformation(:)=0.
        do ipoint=0,Tdomain%n_glob_points-1
            ! recuperation des deformations pour les sorties capteurs
            Tdomain%GrandeurDeformation(ipoint)= Grad_u(0,ipoint) + Grad_u(3,ipoint)
        enddo

    endif

    if (Tdomain%logicD%save_deformation.and.i_snap==0) then
#ifdef MKA3D
        call semname_save_deformation_datdef(nom_dir_sorties,it,Tdomain%mpi_var%my_rank,fnamef)
        open (71,file=fnamef,status="unknown",form="formatted")
        call semname_save_deformation_datepsilon(nom_dir_sorties,it,Tdomain%mpi_var%my_rank,fnamef)
        open (72,file=fnamef,status="unknown",form="formatted")

#else
        call semname_save_deformation_epsilonp(rank,it,fnamef)
        open (71,file=fnamef,status="unknown",form="formatted")
        call semname_save_deformation_epsilonp(rank,it,fnamef)
        open (72,file=fnamef,status="unknown",form="formatted")
#endif


        do ipoint=0,Tdomain%n_glob_points-1
            if (Wheight (ipoint) > 0  ) then
                epsilon_plus  = Grad_u(0,ipoint) + Grad_u(3,ipoint)
                epsilon_minus = (Grad_u(0,ipoint) - Grad_u(3,ipoint))**2 + 4.* Grad_u(1,ipoint) * Grad_u(2,ipoint)
                !!           if ( epsilon_minus < 1.e-60) epsilon_minus = 0. !!initial - gfortran rale
                if ( epsilon_minus < 1.e-20) epsilon_minus = 0.
                epsilon_minus = sqrt (epsilon_minus)

#ifdef MKA3D
                write (71,"(I8,G17.8)") ipoint, epsilon_plus
                write (72,"(I8,G17.8)") ipoint, epsilon_minus
#else
                write (71,"(3G17.8)") Tdomain%GlobCoord (0,ipoint), Tdomain%GlobCoord (1,ipoint), epsilon_plus
                write (72,"(3G17.8)") Tdomain%GlobCoord (0,ipoint), Tdomain%GlobCoord (1,ipoint), epsilon_minus
#endif

            endif

#ifdef MKA3D
            if (Wheight (ipoint) == 0  ) then
                write (71,"(I8,G17.8)") ipoint, zzero
            endif
#endif

        enddo
        call flush(71); call flush (72)
        close(71); close (72)

    endif

    deallocate (Wheight, Grad_u)
    return
end subroutine save_deformation

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
