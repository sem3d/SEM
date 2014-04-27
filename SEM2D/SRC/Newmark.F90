!>
!!\file Newmark.F90
!!\brief Algorithme de Newmark
!!\version 1.0
!!\date 10/03/2009
!! La routine Newmark assure la résolution des équations via un algorithme de predicteur-multi-correcteur
!! des vitesses avec une formulation contrainte-vitesse décalée en temps dans les PML.
!<

module snewmark
    use sdomain
    use scouplage
    use mpi
    implicit none
contains

subroutine Newmark (Tdomain)

    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables
    integer :: ns, ncc,i,j,n,np, ngllx, ngllz, mat, nelem,nf, w_face, nv_aus, nf_aus, nv
    integer :: n_face_pointed, tag_send, tag_receive, i_send, i_stock, ngll, ierr, i_proc
    integer, dimension (MPI_STATUS_SIZE) :: status
    real :: bega, gam1,alpha,dt

    real, dimension (0:1) :: V_free_vertex
    real, dimension (:,:), allocatable :: Vxloc, Vzloc, V_free
    real, dimension (:,:), allocatable :: Smooth


    ! Predictor-MultiCorrector Newmark Velocity Scheme within a
    ! Time staggered Stress-Velocity formulation inside PML

    ! Predictor Phase

    if (Tdomain%TimeD%velocity_scheme) then

        alpha =Tdomain%TimeD%alpha
        bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
        gam1 = 1. / Tdomain%TimeD%gamma

        do n = 0, Tdomain%n_elem -1
            mat = Tdomain%specel(n)%mat_index
            Dt = Tdomain%sSubDomain(mat)%dt

            if (.not. Tdomain%specel(n)%PML) then
                call Prediction_Elem_Veloc (Tdomain%specel(n))
            elseif (Tdomain%specel(n)%CPML) then
                call Prediction_Elem_CPML_Veloc (Tdomain%specel(n),alpha, bega, dt,Vxloc,Vzloc, &
                        Tdomain%sSubDomain(mat)%hPrimez, Tdomain%sSubDomain(mat)%hTPrimex)
            else
                ngllx = Tdomain%specel(n)%ngllx
                ngllz = Tdomain%specel(n)%ngllz
                allocate (Vxloc(0:ngllx-1, 0:ngllz-1))
                allocate (Vzloc(0:ngllx-1, 0:ngllz-1))
                call get_PMLprediction_fv2el (Tdomain,n,Vxloc,vzloc,ngllx,ngllz,alpha, bega,dt)
                if (Tdomain%specel(n)%FPML) then
                    call Prediction_Elem_FPML_Veloc  (Tdomain%specel(n),alpha, bega, dt,Vxloc,Vzloc, &
                        Tdomain%sSubDomain(mat)%hPrimez, Tdomain%sSubDomain(mat)%hTPrimex, &
                        Tdomain%sSubDomain(mat)%freq)
                else
                    call Prediction_Elem_PML_Veloc  (Tdomain%specel(n),alpha, bega, dt,Vxloc,Vzloc, &
                        Tdomain%sSubDomain(mat)%hPrimez, Tdomain%sSubDomain(mat)%hTPrimex)
                endif
                deallocate (VxLoc, Vzloc)
            endif
        enddo

        do n = 0, Tdomain%n_face-1
            mat = Tdomain%sFace(n)%mat_index
            dt = Tdomain%sSubdomain(mat)%dt
            if (.not. Tdomain%sFace(n)%PML)  call Prediction_Face_Veloc (Tdomain%sFace(n))
        enddo

        do n= 0, Tdomain%n_vertex-1
            mat = Tdomain%sVertex(n)%mat_index
            dt = Tdomain%sSubdomain(mat)%dt
            if (.not. Tdomain%sVertex(n)%PML)  call Prediction_Vertex_Veloc (Tdomain%sVertex(n))
        enddo

        ! Solution phase
        do n = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            if (.not. Tdomain%specel(n)%PML ) then
                call get_Displ_fv2el (Tdomain,n)

                call compute_InternalForces_Elem (Tdomain%specel(n),                &
                    Tdomain%sSubDomain(mat)%hprimex,  &
                    Tdomain%sSubDomain(mat)%hTprimex, &
                    Tdomain%sSubDomain(mat)%hprimez,  &
                    Tdomain%sSubDomain(mat)%hTprimez)

            else
                call compute_InternalForces_PML_Elem (Tdomain%specel(n), &
                    Tdomain%sSubDomain(mat)%hprimex, &
                    Tdomain%sSubDomain(mat)%hTprimez)
            endif
        enddo




        ! External Forces
        do n = 0, Tdomain%n_source-1
            do ns =0, Tdomain%sSource(n)%ine-1
                ncc = Tdomain%sSource(n)%Elem(ns)%nr
                ngllx = Tdomain%specel(ncc)%ngllx; ngllz = Tdomain%specel(ncc)%ngllz

                do j = 0,ngllz-1
                    do i = 0,ngllx-1
                        do np = 0,1
                            Tdomain%specel(ncc)%Forces(i,j,np) = Tdomain%specel(ncc)%Forces(i,j,np) +   &
                                CompSource (Tdomain%sSource(n),Tdomain%TimeD%rtime)*  Tdomain%sSource(n)%Elem(ns)%ExtForce(i,j,np)
                        enddo
                    enddo
                enddo

            enddo
        enddo

        ! Communication of Forces

        do nf = 0, Tdomain%n_face-1
            Tdomain%sFace(nf)%Forces = 0
            if  (Tdomain%sFace(nf)%PML) then
                Tdomain%sFace(nf)%Forces1 = 0; Tdomain%sFace(nf)%Forces2 = 0
            endif
            nelem = Tdomain%sFace(nf)%Near_element(0)
            w_face = Tdomain%sFace(nf)%Which_face(0)
            call getInternalF_el2f (Tdomain,nelem,nf,w_face,.true.)
            if (Tdomain%sFace(nf)%PML ) call getInternalF_PML_el2f (Tdomain,nelem,nf,w_face,.true.)
            nelem = Tdomain%sFace(nf)%Near_element(1)
            if (nelem > -1) then
                w_face = Tdomain%sFace(nf)%Which_face(1)
                call getInternalF_el2f (Tdomain,nelem,nf,w_face,Tdomain%sFace(nf)%coherency)
                if (Tdomain%sFace(nf)%PML ) call getInternalF_PML_el2f (Tdomain,nelem,nf,w_face,Tdomain%sFace(nf)%coherency)
            endif
        enddo

        do nv = 0, Tdomain%n_vertex-1
            Tdomain%sVertex(nv)%Forces= 0
            if (Tdomain%sVertex(nv)%PML ) then
                Tdomain%sVertex(nv)%Forces1= 0;    Tdomain%sVertex(nv)%Forces2= 0
            endif
        enddo



        do n = 0, Tdomain%n_elem-1
            ngllx= Tdomain%specel(n)%ngllx;     ngllz= Tdomain%specel(n)%ngllz


            nv = Tdomain%specel(n)%Near_Vertex(0)
            Tdomain%sVertex(nv)%Forces(0:1) = Tdomain%sVertex(nv)%Forces(0:1) + Tdomain%specel(n)%Forces(0,0,0:1)

            if (Tdomain%sVertex(nv)%PML) then
                Tdomain%sVertex(nv)%Forces1(0:1) = Tdomain%sVertex(nv)%Forces1(0:1) + Tdomain%specel(n)%Forces1(0,0,0:1)
                Tdomain%sVertex(nv)%Forces2(0:1) = Tdomain%sVertex(nv)%Forces2(0:1) + Tdomain%specel(n)%Forces2(0,0,0:1)
            endif


            nv = Tdomain%specel(n)%Near_Vertex(1)
            Tdomain%sVertex(nv)%Forces(0:1) = Tdomain%sVertex(nv)%Forces(0:1) + Tdomain%specel(n)%Forces(ngllx-1,0,0:1)
            if (Tdomain%sVertex(nv)%PML) then
                Tdomain%sVertex(nv)%Forces1(0:1) = Tdomain%sVertex(nv)%Forces1(0:1) + Tdomain%specel(n)%Forces1(ngllx-1,0,0:1)
                Tdomain%sVertex(nv)%Forces2(0:1) = Tdomain%sVertex(nv)%Forces2(0:1) + Tdomain%specel(n)%Forces2(ngllx-1,0,0:1)
            endif

            nv = Tdomain%specel(n)%Near_Vertex(2)
            Tdomain%sVertex(nv)%Forces(0:1) = Tdomain%sVertex(nv)%Forces(0:1) + Tdomain%specel(n)%Forces(ngllx-1,ngllz-1,0:1)
            if (Tdomain%sVertex(nv)%PML) then
                Tdomain%sVertex(nv)%Forces1(0:1) = Tdomain%sVertex(nv)%Forces1(0:1) + Tdomain%specel(n)%Forces1(ngllx-1,ngllz-1,0:1)
                Tdomain%sVertex(nv)%Forces2(0:1) = Tdomain%sVertex(nv)%Forces2(0:1) + Tdomain%specel(n)%Forces2(ngllx-1,ngllz-1,0:1)
            endif

            nv = Tdomain%specel(n)%Near_Vertex(3)
            Tdomain%sVertex(nv)%Forces(0:1) = Tdomain%sVertex(nv)%Forces(0:1) + Tdomain%specel(n)%Forces(0,ngllz-1,0:1)
            if (Tdomain%sVertex(nv)%PML) then
                Tdomain%sVertex(nv)%Forces1(0:1) = Tdomain%sVertex(nv)%Forces1(0:1) + Tdomain%specel(n)%Forces1(0,ngllz-1,0:1)
                Tdomain%sVertex(nv)%Forces2(0:1) = Tdomain%sVertex(nv)%Forces2(0:1) + Tdomain%specel(n)%Forces2(0,ngllz-1,0:1)
            endif


        enddo




        ! AJOUT DES FORCES MKA3D
#ifdef COUPLAGE
        if (Tdomain%TimeD%ntime>0) then
            call calcul_couplage_force(Tdomain,Tdomain%TimeD%ntime)
        endif

        ! la ForcesMka corespond a la contrainte multiplie par la surface du point de gauss correspondant
        ! sur un meme proc la somme est deja faite
        ! par contre lorsqu un vertex est partage sur plusieurs proc alors chaque proc n a qu une partie de la somme
        ! il faut donc lui ajouter les contributions des autres proc
        ! pour prendre en compte les forces imposee lors du couplage avec mka sur les points de gauss internes aux faces
        do nf = 0, Tdomain%n_face-1
            ngll = Tdomain%sFace(nf)%ngll
            do i=1,ngll-2
                Tdomain%sFace(nf)%Forces(i,0:1) = Tdomain%sFace(nf)%ForcesMka(i,0:1) + Tdomain%sFace(nf)%Forces(i,0:1)
            enddo

        enddo


        ! pour prendre en compte les forces imposee lors du couplage avec mka sur les points de gauss des vertex
        do nv = 0, Tdomain%n_vertex-1
            Tdomain%sVertex(nv)%Forces(0:1) = Tdomain%sVertex(nv)%ForcesMka(0:1) + Tdomain%sVertex(nv)%Forces(0:1)
        enddo

#endif





        ! Communicate forces among processors

        ! Double value on the vertices
        do n = 0, Tdomain%n_communications - 1
            do nv = 0, Tdomain%sWall(n)%n_vertices-1
                nv_aus = Tdomain%sWall(n)%Vertex_List(nv)
                Tdomain%sVertex(nv_aus)%Double_Value(0:1) = Tdomain%sVertex(nv_aus)%Forces(0:1)
            enddo
        enddo


        !==========================================================
        do i_proc = 0, Tdomain%n_communications - 1
            i_send = Tdomain%Communication_list (i_proc)
            i_stock = 0

            do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
                n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
                ngll = Tdomain%sFace(n_face_pointed)%ngll




                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    Tdomain%sWall(i_proc)%Send_data_2(i_stock:i_stock+ngll-3,0:1) = Tdomain%sFace(n_face_pointed)%Forces (1:ngll-2,0:1)
                else
                    do j = 1, ngll-2
                        Tdomain%sWall(i_proc)%Send_data_2(i_stock+j-1,0:1) = Tdomain%sFace(n_face_pointed)%Forces ( ngll-1-j,0:1)
                    enddo
                endif
                i_stock = i_stock + ngll - 2
            enddo



            do nv = 0, Tdomain%sWall(i_proc)%n_vertices - 1
                nv_aus =  Tdomain%sWall(i_proc)%Vertex_List(nv)
                Tdomain%sWall(i_proc)%Send_data_2 (i_stock,0:1)  = Tdomain%sVertex(nv_aus)%Double_Value(0:1)
                i_stock = i_stock + 1
            enddo


            do nf = 0, Tdomain%sWall(i_proc)%n_pml_faces - 1
                n_face_pointed = Tdomain%sWall(i_proc)%FacePML_List(nf)
                ngll = Tdomain%sFace(n_face_pointed)%ngll
                if (Tdomain%sWall(i_proc)%FacePML_Coherency(nf)) then
                    Tdomain%sWall(i_proc)%Send_data_2(i_stock:i_stock+ngll-3,0:1) = Tdomain%sFace(n_face_pointed)%Forces1 (1:ngll-2,0:1)
                    i_stock = i_stock + ngll - 2
                    Tdomain%sWall(i_proc)%Send_data_2(i_stock:i_stock+ngll-3,0:1) = Tdomain%sFace(n_face_pointed)%Forces2 (1:ngll-2,0:1)
                else
                    do j = 1, ngll-2
                        Tdomain%sWall(i_proc)%Send_data_2(i_stock+j-1,0:1) = Tdomain%sFace(n_face_pointed)%Forces1 ( ngll-1-j,0:1)
                    enddo
                    i_stock = i_stock + ngll - 2
                    do j = 1, ngll-2
                        Tdomain%sWall(i_proc)%Send_data_2(i_stock+j-1,0:1) = Tdomain%sFace(n_face_pointed)%Forces2 ( ngll-1-j,0:1)
                    enddo
                endif
                i_stock = i_stock + ngll - 2
            enddo

            tag_send = i_send * Tdomain%MPI_var%n_proc +Tdomain%MPI_var%my_rank + 900
            tag_receive = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 900

            call MPI_SEND (Tdomain%sWall(i_proc)%Send_data_2,2* Tdomain%sWall(i_proc)%n_points, MPI_DOUBLE_PRECISION, i_send, &
                tag_send, Tdomain%communicateur, ierr )

            call MPI_RECV (Tdomain%sWall(i_proc)%Receive_data_2, 2* Tdomain%sWall(i_proc)%n_points, MPI_DOUBLE_PRECISION, i_send, &
                tag_receive, Tdomain%communicateur, status, ierr )

            i_stock = 0
            do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
                n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
                ngll = Tdomain%sFace(n_face_pointed)%ngll


                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    Tdomain%sFace(n_face_pointed)%Forces (1:ngll-2,0:1) =  Tdomain%sFace(n_face_pointed)%Forces ( 1:ngll-2,0:1)  + &
                        Tdomain%sWall(i_proc)%Receive_data_2(i_stock:i_stock+ngll-3,0:1)
                else
                    do j = 1, ngll-2
                        Tdomain%sFace(n_face_pointed)%Forces (ngll-1-j,0:1) =  Tdomain%sFace(n_face_pointed)%Forces (ngll-1-j,0:1)+ &
                            Tdomain%sWall(i_proc)%Receive_data_2(i_stock+j-1,0:1)

                    enddo
                endif
                i_stock = i_stock + ngll - 2
            enddo


            do nv = 0, Tdomain%sWall(i_proc)%n_vertices - 1
                nv_aus =  Tdomain%sWall(i_proc)%Vertex_List(nv)
                Tdomain%sVertex(nv_aus)%Forces(0:1) = Tdomain%sVertex(nv_aus)%Forces(0:1) + &
                    Tdomain%sWall(i_proc)%Receive_data_2 (i_stock,0:1)
                i_stock = i_stock + 1
            enddo






            do nf = 0, Tdomain%sWall(i_proc)%n_pml_faces - 1
                n_face_pointed = Tdomain%sWall(i_proc)%FacePML_List(nf)
                ngll = Tdomain%sFace(n_face_pointed)%ngll
                if (Tdomain%sWall(i_proc)%FacePML_Coherency(nf)) then
                    Tdomain%sFace(n_face_pointed)%Forces1 (1:ngll-2,0:1) =  Tdomain%sFace(n_face_pointed)%Forces1 ( 1:ngll-2,0:1)  + &
                        Tdomain%sWall(i_proc)%Receive_data_2(i_stock:i_stock+ngll-3,0:1)
                    i_stock = i_stock + ngll - 2
                    Tdomain%sFace(n_face_pointed)%Forces2 (1:ngll-2,0:1) =  Tdomain%sFace(n_face_pointed)%Forces2 ( 1:ngll-2,0:1)  + &
                        Tdomain%sWall(i_proc)%Receive_data_2(i_stock:i_stock+ngll-3,0:1)
                else
                    do j = 1, ngll-2
                        Tdomain%sFace(n_face_pointed)%Forces1 (ngll-1-j,0:1) =  Tdomain%sFace(n_face_pointed)%Forces1 (ngll-1-j,0:1)+ &
                            Tdomain%sWall(i_proc)%Receive_data_2(i_stock+j-1,0:1)
                    enddo
                    i_stock = i_stock + ngll - 2
                    do j = 1, ngll-2
                        Tdomain%sFace(n_face_pointed)%Forces2 (ngll-1-j,0:1) =  Tdomain%sFace(n_face_pointed)%Forces2 (ngll-1-j,0:1)+ &
                            Tdomain%sWall(i_proc)%Receive_data_2(i_stock+j-1,0:1)
                    enddo
                endif
                i_stock = i_stock + ngll - 2
            enddo
        enddo

        !========fin des comm entre proc=============


        !==========================================================
        ! Definition of a super object

        if (Tdomain%logicD%super_object_local_present) then
            do n = 0, Tdomain%n_fault-1
                do nf = 0, Tdomain%sFault(n)%n_face-1
                    ngllx = Tdomain%sFault(n)%fFace(nf)%ngll
                    mat = Tdomain%sFault(n)%fFace(nf)%mat_index
                    allocate (V_free(1:ngllx-2,0:1))
                    nf_aus = Tdomain%sFault(n)%fFace(nf)%Face_UP
                    call get_vfree_face(Tdomain%sFace(nf_aus),V_free,ngllx,.false.,.true.)
                    nf_aus = Tdomain%sFault(n)%fFace(nf)%Face_DOWN
                    call get_vfree_face(Tdomain%sFace(nf_aus),V_free,ngllx,.true.,Tdomain%sFault(n)%fFace(nf)%Coherency)
                    select case (Tdomain%sFault(n)%Problem_type)
                    case (1)
                        call traction_on_face_sw (Tdomain%sFault(n)%fFace(nf), V_free,Tdomain%sSubdomain(mat)%dt, &
                            Tdomain%sFault(n)%imposed_Tolerance)
                    case (2)
                        call traction_on_face_adhesion (Tdomain%sFault(n)%fFace(nf), V_free,Tdomain%sSubdomain(mat)%dt, &
                            Tdomain%sFault(n)%imposed_Tolerance)
                    end select
                    deallocate (V_free)
                enddo

                do nf = 0, Tdomain%sFault(n)%n_vertex-1
                    if (.not. Tdomain%sFault(n)%fVertex(nf)%Termination) then
                        mat = Tdomain%sFault(n)%fVertex(nf)%mat_index
                        nv_aus = Tdomain%sFault(n)%fVertex(nf)%Vertex_UP
                        call get_vfree_vertex (Tdomain%sVertex(nv_aus), V_free_Vertex, .false.)
                        nv_aus = Tdomain%sFault(n)%fVertex(nf)%Vertex_DOWN
                        call get_vfree_vertex (Tdomain%sVertex(nv_aus), V_free_Vertex, .true.)
                        select case (Tdomain%sFault(n)%Problem_type)
                        case (1)
                            call traction_on_vertex_sw (Tdomain%sFault(n)%fVertex(nf),    &
                                V_free_vertex,Tdomain%sSubdomain(mat)%dt, Tdomain%sFault(n)%imposed_Tolerance)

                        case (2)
                            call traction_on_vertex_adhesion (Tdomain%sFault(n)%fVertex(nf),   &
                                V_free_vertex,Tdomain%sSubdomain(mat)%dt, Tdomain%sFault(n)%imposed_Tolerance)
                        end select
                    endif
                enddo
            enddo




            ! Smoothing procedure for the normal traction
            do  n = 0, Tdomain%n_fault-1
                if (Tdomain%sFault(n)%smoothing ) then
                    ns = -1
                    do nf = 0, Tdomain%sFault(n)%n_face-1
                        ns = ns + Tdomain%sFault(n)%fFace(nf)%ngll-1
                    enddo
                    allocate (Smooth(0:ns-1,0:1))
                    ncc = 0
                    do nf = 0, Tdomain%sFault(n)%n_face-1
                        ngllx = Tdomain%sFault(n)%fFace(nf)%ngll
                        Smooth (ncc:ncc+ngllx-3,0) = Tdomain%sFault(n)%fFace(nf)%distance(1:ngllx-2)
                        Smooth (ncc:ncc+ngllx-3,1) = Tdomain%sFault(n)%fFace(nf)%traction(1:ngllx-2,1)
                        nv_aus = Tdomain%sFault(n)%fFace(nf)%Face_to_Vertex(1)
                        if (.not. Tdomain%sFault(n)%fVertex(nv_aus)%Termination)  then
                            Smooth(ncc+ngllx-2,0) = Tdomain%sFault(n)%fVertex(nv_aus)%distance
                            Smooth(ncc+ngllx-2,1) = Tdomain%sFault(n)%fVertex(nv_aus)%traction(1)
                        endif
                        ncc = ncc + ngllx- 1
                    enddo

                    !       call smoothing_exp(smooth,ns-1,1,4,Tdomain%sFault(n)%dx_smoothing)
                    call smoothing_exp(smooth,ns-1,1,4,0.4)
                    ncc = 0
                    do nf = 0, Tdomain%sFault(n)%n_face-1
                        ngllx = Tdomain%sFault(n)%fFace(nf)%ngll
                        Tdomain%sFault(n)%fFace(nf)%traction(1:ngllx-2,1) = Smooth (ncc:ncc+ngllx-3,1)
                        nv_aus = Tdomain%sFault(n)%fFace(nf)%Face_to_Vertex(1)
                        if (.not. Tdomain%sFault(n)%fVertex(nv_aus)%Termination)  &
                            Tdomain%sFault(n)%fVertex(nv_aus)%traction(1) = Smooth(ncc+ngllx-2,1)
                        ncc = ncc + ngllx- 1
                    enddo
                    deallocate (Smooth)
                    print *, "here you are doing smmothing"
                endif

                do nf = 0, Tdomain%sFault(n)%n_face-1
                    call rotate_traction_face (Tdomain%sFault(n)%fFace(nf))
                enddo
                do nf = 0, Tdomain%sFault(n)%n_vertex-1
                    call rotate_traction_vertex (Tdomain%sFault(n)%fVertex(nf))
                enddo
            enddo


            ! Following the smoothing - special treatement

            do n = 0, Tdomain%n_fault-1
                do nf = 0, Tdomain%sFault(n)%n_face-1
                    ngllx = Tdomain%sFault(n)%fFace(nf)%ngll
                    nf_aus = Tdomain%sFault(n)%fFace(nf)%Face_UP
                    do i = 0,1
                        Tdomain%sFace(nf_aus)%Forces(1:ngllx-2,i) = Tdomain%sFace(nf_aus)%Forces(1:ngllx-2,i) + &
                            Tdomain%sFault(n)%fFace(nf)%Bt(1:ngllx-2) * Tdomain%sFault(n)%fFace(nf)%Traction(1:ngllx-2,i)
                    enddo
                    nf_aus = Tdomain%sFault(n)%fFace(nf)%Face_DOWN
                    if (Tdomain%sFault(n)%fFace(nf)%Coherency) then
                        do i =0,1
                            Tdomain%sFace(nf_aus)%Forces(1:ngllx-2,i) = Tdomain%sFace(nf_aus)%Forces(1:ngllx-2,i) - &
                                Tdomain%sFault(n)%fFace(nf)%Bt(1:ngllx-2) * Tdomain%sFault(n)%fFace(nf)%Traction(1:ngllx-2,i)
                        enddo
                    else
                        do i =0,1
                            do j = 1,ngllx-2
                                Tdomain%sFace(nf_aus)%Forces(j,i) = Tdomain%sFace(nf_aus)%Forces(j,i) - &
                                    Tdomain%sFault(n)%fFace(nf)%Bt(ngllx-1-j) * Tdomain%sFault(n)%fFace(nf)%Traction(ngllx-1-j,i)
                            enddo
                        enddo
                    endif
                enddo

                do nf = 0, Tdomain%sFault(n)%n_vertex-1
                    if (.not. Tdomain%sFault(n)%fVertex(nf)%Termination) then
                        nv_aus = Tdomain%sFault(n)%fVertex(nf)%Vertex_UP
                        Tdomain%sVertex(nv_aus)%Forces(0:1) = Tdomain%sVertex(nv_aus)%Forces(0:1) + &
                            Tdomain%sFault(n)%fVertex(nf)%Bt * Tdomain%sFault(n)%fVertex(nf)%Traction(0:1)
                        nv_aus = Tdomain%sFault(n)%fVertex(nf)%Vertex_DOWN
                        Tdomain%sVertex(nv_aus)%Forces(0:1) = Tdomain%sVertex(nv_aus)%Forces(0:1) - &
                            Tdomain%sFault(n)%fVertex(nf)%Bt * Tdomain%sFault(n)%fVertex(nf)%Traction(0:1)
                    endif
                enddo
            enddo

        endif

        !===fin traitement super objet =============


        do n = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            if (.not. Tdomain%specel(n)%PML) then
                call Correction_Elem_Veloc (Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt)
            else
                if (Tdomain%specel(n)%FPML) then
                    call Correction_Elem_FPML_Veloc (Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt, Tdomain%sSubdomain(mat)%freq)
                else
                    call Correction_Elem_PML_Veloc (Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt)
                endif
            endif
        enddo

        do nf = 0, Tdomain%n_face-1
            mat = Tdomain%sFace(nf)%mat_index
            if (.not. Tdomain%sFace(nf)%PML) then
                call Correction_Face_Veloc (Tdomain%sFace(nf),Tdomain%sFace(nf)%ngll, Tdomain%sSubDomain(mat)%Dt)
            else
                if (Tdomain%sFace(nf)%FPML) then
                    call Correction_Face_FPML_Veloc (Tdomain%sFace(nf),Tdomain%sSubDomain(mat)%Dt, Tdomain%sSubDomain(mat)%freq)
                else
                    call Correction_Face_PML_Veloc (Tdomain%sFace(nf),Tdomain%sSubDomain(mat)%Dt)
                endif
            endif
        enddo

        do nv= 0, Tdomain%n_vertex-1
            mat = Tdomain%sVertex(nv)%mat_index
            if (.not. Tdomain%sVertex(nv)%PML) then
                call Correction_Vertex_Veloc (Tdomain%sVertex(nv), Tdomain%sSubDomain(mat)%Dt)
            else
                if (Tdomain%sVertex(nv)%FPML) then
                    call Correction_Vertex_FPML_Veloc (Tdomain%sVertex(nv),Tdomain%sSubDomain(mat)%Dt, Tdomain%sSubdomain(mat)%freq)
                else
                    call Correction_Vertex_PML_Veloc (Tdomain%sVertex(nv),Tdomain%sSubDomain(mat)%Dt)
                endif
            endif
        enddo

    endif   ! if Velocity Scheme




    return
end subroutine Newmark

end module snewmark
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
