!>
!!\file define_fault_properties.F90
!!\brief Definition des defauts
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type(Domain), intent (INOUT) Tdomain
!<


subroutine define_fault_properties (Tdomain)

    ! Modified by Gaetano Festa 01/06/2005

    use sdomain
    use semdatafiles
    use mpi

    implicit none
    type(Domain), intent (INOUT) :: Tdomain

    integer :: n,nf, n_begin, n_end, nv, j, ngll, n_anomalies, k, i, i_proc_of_anomaly
    integer :: i_send, tag_send, tag_receive, ierr
    integer, dimension (MPI_STATUS_SIZE) :: status

    real :: sigma1, sigma3, angle, mus, mud, Dc, pig, tau0_percentage
    real :: w,b,Ct, Deltau_res,c1,s1, sxx, sxy, syy
    real, dimension (:), allocatable :: Send_data, Receive_data

    character (len=MAX_FILE_SIZE) :: fnamef

    nf = 0
    do n = 0, Tdomain%n_super_object-1
        if (Tdomain%Super_object_type(n) == "F" ) then

            ! Communicating normal between points
            do j = 0, Tdomain%n_communications -1
                i_send = Tdomain%Communication_List(j)

                if (Tdomain%sWall(j)%n_vertex_superobject > 0) then
                    allocate (Send_Data(0:Tdomain%sWall(j)%n_vertex_superobject*2-1))
                    allocate (Receive_Data(0:Tdomain%sWall(j)%n_vertex_superobject*2-1))
                    do i = 0, Tdomain%sWall(j)%n_vertex_superobject-1
                        k = Tdomain%sWall(j)%Vertex_SuperObject_List(i)
                        Send_data(2*i:2*i+1) = Tdomain%sFault(nf)%fVertex(k)%normal(0:1)
                    enddo

                    tag_send = i_send * Tdomain%MPI_var%n_proc +Tdomain%MPI_var%my_rank + 800
                    tag_receive = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 800

                    call MPI_SEND (Send_data,2*Tdomain%sWall(j)%n_vertex_superobject, MPI_DOUBLE_PRECISION, i_send, &
                        tag_send, Tdomain%communicateur, ierr )
                    call MPI_RECV (Receive_data, 2*Tdomain%sWall(j)%n_vertex_superobject, MPI_DOUBLE_PRECISION, i_send, &
                        tag_receive, Tdomain%communicateur, status, ierr )

                    do i = 0, Tdomain%sWall(j)%n_vertex_superobject-1
                        k = Tdomain%sWall(j)%Vertex_SuperObject_List(i)
                        Tdomain%sFault(nf)%fVertex(k)%normal(0) = 0.5 * (Tdomain%sFault(nf)%fVertex(k)%normal(0) + Receive_data(2*i))
                        Tdomain%sFault(nf)%fVertex(k)%normal(1) = 0.5 * (Tdomain%sFault(nf)%fVertex(k)%normal(1) + Receive_data(2*i+1))
                    enddo
                    deallocate (Send_data, Receive_data)
                endif
            enddo

            call semname_define_fault_data(Tdomain%super_object_file(n),fnamef)
            open (24,file=fnamef,status="old", form="formatted")

            read (24,*)
            read (24,*) Tdomain%sFault(nf)%Problem_type
            read (24,*) Tdomain%sFault(nf)%imposed_tolerance
            read (24,*) Tdomain%logicD%save_fault_trace
            read (24,*) Tdomain%sFault(nf)%smoothing
            if (Tdomain%sFault(nf)%smoothing) then
                read (24,*) Tdomain%sFault(nf)%dx_smoothing
                write (*,*) "Here you choose a smoothing of the normal traction"
                write (*,*) "points on the fault are required to be ordered"
            else
                read(24,*)
            endif

            if (Tdomain%sFault(nf)%Problem_type == 1 ) then

                read (24,*) sigma1
                read (24,*) sigma3
                read (24,*) angle
                read (24,*) mus
                read (24,*) mud
                read (24,*) Dc
                pig = acos (-1.)
                angle = angle *pig/180.

                ! Transform orignal stress in carthesian coordinates
                c1 = cos(angle)
                s1 = sin(angle)
                sxx = c1**2 * sigma1 + s1**2 * sigma3
                syy = c1**2 * sigma3 + s1**2 * sigma1
                sxy = c1*s1*(sigma3-sigma1)
                do j = 0, Tdomain%sFault(nf)%n_vertex-1
                    Tdomain%sFault(nf)%fVertex(j)%tau0 =  sxx * Tdomain%sFault(nf)%fVertex(j)%normal(0) &
                        + sxy * Tdomain%sFault(nf)%fVertex(j)%normal(1)
                    Tdomain%sFault(nf)%fVertex(j)%sigma0 =  sxy * Tdomain%sFault(nf)%fVertex(j)%normal(0) &
                        + syy * Tdomain%sFault(nf)%fVertex(j)%normal(1)
                    Tdomain%sFault(nf)%fVertex(j)%mus = mus
                    Tdomain%sFault(nf)%fVertex(j)%mud = mud
                    Tdomain%sFault(nf)%fVertex(j)%Dc = Dc
                    Tdomain%sFault(nf)%fVertex(j)%Deltav = 0
                    Tdomain%sFault(nf)%fVertex(j)%Deltau = 0
                    Tdomain%sFault(nf)%fVertex(j)%AbsDeltau0 = 0.
                enddo
                do j = 0, Tdomain%sFault(nf)%n_face-1
                    ngll = Tdomain%sFault(nf)%fFace(j)%ngll
                    allocate (Tdomain%sFault(nf)%fFace(j)%tau0 (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%sigma0 (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%mus (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%mud (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Dc (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Bt(0:ngll-1))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Cgamma(1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Deltav(1:ngll-2,0:1))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Deltau(1:ngll-2,0:1))
                    allocate (Tdomain%sFault(nf)%fFace(j)%AbsDeltau0(1:ngll-2))

                    Tdomain%sFault(nf)%fFace(j)%Deltav = 0
                    Tdomain%sFault(nf)%fFace(j)%Deltau = 0
                    Tdomain%sFault(nf)%fFace(j)%AbsDeltau0 = 0.
                    allocate (Tdomain%sFault(nf)%fFace(j)%Traction(1:ngll-2,0:1))
                    do  i = 1, ngll-2
                        Tdomain%sFault(nf)%fFace(j)%sigma0(i) = sxy * Tdomain%sFault(nf)%fFace(j)%normal(i,0) &
                            + syy * Tdomain%sFault(nf)%fFace(j)%normal(i,1)
                        Tdomain%sFault(nf)%fFace(j)%tau0(i) =  sxx * Tdomain%sFault(nf)%fFace(j)%normal(i,0) &
                            + sxy * Tdomain%sFault(nf)%fFace(j)%normal(i,1)
                        Tdomain%sFault(nf)%fFace(j)%mus(i) = mus
                        Tdomain%sFault(nf)%fFace(j)%mud(i) = mud
                        Tdomain%sFault(nf)%fFace(j)%Dc(i) = Dc
                    enddo
                enddo

                ! Rotate along the normal - tangential reference frame
                ! (t,n) is a leftward reference t = (-n2, n1) ; n = (n1, n2)

                do j = 0, Tdomain%sFault(nf)%n_vertex-1
                    c1 = - Tdomain%sFault(nf)%fVertex(j)%tau0 * Tdomain%sFault(nf)%fVertex(j)%normal(1) + &
                        Tdomain%sFault(nf)%fVertex(j)%sigma0 * Tdomain%sFault(nf)%fVertex(j)%normal(0)
                    s1 = Tdomain%sFault(nf)%fVertex(j)%tau0 * Tdomain%sFault(nf)%fVertex(j)%normal(0) + &
                        Tdomain%sFault(nf)%fVertex(j)%sigma0 * Tdomain%sFault(nf)%fVertex(j)%normal(1)
                    Tdomain%sFault(nf)%fVertex(j)%tau0 = c1
                    Tdomain%sFault(nf)%fVertex(j)%sigma0 =  s1
                enddo
                do j = 0, Tdomain%sFault(nf)%n_face-1
                    ngll = Tdomain%sFault(nf)%fFace(j)%ngll
                    do  i = 1, ngll-2
                        c1 = - Tdomain%sFault(nf)%fFace(j)%tau0(i) * Tdomain%sFault(nf)%fFace(j)%normal(i,1) + &
                            Tdomain%sFault(nf)%fFace(j)%sigma0(i) * Tdomain%sFault(nf)%fFace(j)%normal(i,0)
                        s1 = Tdomain%sFault(nf)%fFace(j)%tau0(i) * Tdomain%sFault(nf)%fFace(j)%normal(i,0) + &
                            Tdomain%sFault(nf)%fFace(j)%sigma0(i) * Tdomain%sFault(nf)%fFace(j)%normal(i,1)
                        Tdomain%sFault(nf)%fFace(j)%tau0(i) = c1
                        Tdomain%sFault(nf)%fFace(j)%sigma0(i) = s1
                    enddo
                enddo
                read (24,*) n_anomalies
                do k = 0, n_anomalies-1
                    read (24,*) i_proc_of_anomaly

                    if (i_proc_of_anomaly == Tdomain%MPI_var%my_rank) then
                        read (24,*) n_begin
                        read (24,*) n_end
                        read (24,*) tau0_percentage
                        read (24,*) mus
                        read (24,*) mud
                        read (24,*) Dc
                        tau0_percentage = tau0_percentage/100.

                        do j = n_begin, n_end
                            Tdomain%sFault(nf)%fFace(j)%tau0 = Tdomain%sFault(nf)%fFace(j)%tau0 * (1.+tau0_percentage)* &
                                abs(Tdomain%sFault(nf)%fFace(j)%sigma0/Tdomain%sFault(nf)%fFace(j)%tau0) * mus
                            Tdomain%sFault(nf)%fFace(j)%mus = mus
                            Tdomain%sFault(nf)%fFace(j)%mud = mud
                            Tdomain%sFault(nf)%fFace(j)%Dc = Dc
                            nv = Tdomain%sFault(nf)%fFace(j)%Face_to_Vertex(0)
                            Tdomain%sFault(nf)%fVertex(nv)%tau0 =  Tdomain%sFault(nf)%fVertex(nv)%tau0 * (1.+tau0_percentage)* &
                                abs(Tdomain%sFault(nf)%fVertex(nv)%sigma0/Tdomain%sFault(nf)%fVertex(nv)%tau0) * mus
                            Tdomain%sFault(nf)%fVertex(nv)%mus = mus
                            Tdomain%sFault(nf)%fVertex(nv)%mud = mud
                            Tdomain%sFault(nf)%fVertex(nv)%Dc = Dc
                            nv = Tdomain%sFault(nf)%fFace(j)%Face_to_Vertex(1)
                            Tdomain%sFault(nf)%fVertex(nv)%tau0 =  Tdomain%sFault(nf)%fVertex(nv)%tau0 * (1.+tau0_percentage)* &
                                abs(Tdomain%sFault(nf)%fVertex(nv)%sigma0/Tdomain%sFault(nf)%fVertex(nv)%tau0) * mus
                            Tdomain%sFault(nf)%fVertex(nv)%mus = mus
                            Tdomain%sFault(nf)%fVertex(nv)%mud = mud
                            Tdomain%sFault(nf)%fVertex(nv)%Dc = Dc
                        enddo
                    else
                        do j = 1, 6
                            read (24,*)
                        enddo
                    endif
                enddo
                close (24)

            else if (Tdomain%sFault(nf)%Problem_type == 2 ) then
                read (24,*) sigma1
                read (24,*) sigma3
                read (24,*) angle
                read (24,*) mud
                read (24,*) Ct
                read (24,*) w
                read (24,*) b
                read (24,*) Deltau_res
                do j = 0, Tdomain%sFault(nf)%n_vertex-1
                    c1 = cos(angle)*Tdomain%sFault(nf)%fVertex(j)%normal(1) -sin(angle)*Tdomain%sFault(nf)%fVertex(j)%normal(0)
                    s1 = sin(angle)*Tdomain%sFault(nf)%fVertex(j)%normal(1) +cos(angle)*Tdomain%sFault(nf)%fVertex(j)%normal(0)
                    Tdomain%sFault(nf)%fVertex(j)%tau0 = c1*s1*(sigma3-sigma1)
                    Tdomain%sFault(nf)%fVertex(j)%sigma0 = c1**2*sigma3+s1**2*sigma1
                    Tdomain%sFault(nf)%fVertex(j)%mud = mud
                    Tdomain%sFault(nf)%fVertex(j)%w = w
                    Tdomain%sFault(nf)%fVertex(j)%Ct = Ct
                    Tdomain%sFault(nf)%fVertex(j)%b = b
                    Tdomain%sFault(nf)%fVertex(j)%beta_p = 1
                    Tdomain%sFault(nf)%fVertex(j)%log_vertex = .true.
                    Tdomain%sFault(nf)%fVertex(j)%Deltav = 0
                    Tdomain%sFault(nf)%fVertex(j)%Deltau(0) = -Deltau_res
                    Tdomain%sFault(nf)%fVertex(j)%Deltau(1) = 0
                enddo
                do j = 0, Tdomain%sFault(nf)%n_face-1
                    ngll = Tdomain%sFault(nf)%fFace(j)%ngll
                    allocate (Tdomain%sFault(nf)%fFace(j)%tau0 (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%sigma0 (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Ct (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%b (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%w (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%beta_p (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%mud (1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Bt(0:ngll-1))
                    allocate (Tdomain%sFault(nf)%fFace(j)%log_face(1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Cgamma(1:ngll-2))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Deltav(1:ngll-2,0:1))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Deltau(1:ngll-2,0:1))
                    allocate (Tdomain%sFault(nf)%fFace(j)%Traction(1:ngll-2,0:1))

                    do  i = 1, ngll-2
                        c1 = cos(angle)*Tdomain%sFault(nf)%fFace(j)%normal(i,1) -sin(angle)*Tdomain%sFault(nf)%fFace(j)%normal(i,0)
                        s1 = sin(angle)*Tdomain%sFault(nf)%fFace(j)%normal(i,1) +cos(angle)*Tdomain%sFault(nf)%fFace(j)%normal(i,0)
                        Tdomain%sFault(nf)%fFace(j)%sigma0(i) = c1**2*sigma3+s1**2*sigma1
                        Tdomain%sFault(nf)%fFace(j)%tau0(i) = c1*s1*(sigma3-sigma1)
                        Tdomain%sFault(nf)%fFace(j)%mud(i) = mud
                        Tdomain%sFault(nf)%fFace(j)%Ct(i) = Ct
                        Tdomain%sFault(nf)%fFace(j)%b(i) = b
                        Tdomain%sFault(nf)%fFace(j)%w(i) = w
                        Tdomain%sFault(nf)%fFace(j)%beta_p(i) = 1.
                        Tdomain%sFault(nf)%fFace(j)%log_face(i) = .true.
                        Tdomain%sFault(nf)%fFace(j)%Deltav(i,:) = 0
                        Tdomain%sFault(nf)%fFace(j)%Deltau(i,0) = - Deltau_res
                        Tdomain%sFault(nf)%fFace(j)%Deltau(i,1) = 0
                    enddo
                enddo
                read (24,*) n_anomalies
                do k = 0, n_anomalies-1
                    read (24,*) i_proc_of_anomaly

                    if (i_proc_of_anomaly == Tdomain%MPI_var%my_rank) then
                        read (24,*) n_begin
                        read (24,*) n_end
                        read (24,*) tau0_percentage
                        read (24,*) mud
                        read (24,*) Ct
                        read (24,*) w
                        read (24,*) b
                        read (24,*) Deltau_res
                        do j = n_begin, n_end

                            if (c1*s1*(sigma3-sigma1) > 0.) then
                                Tdomain%sFault(nf)%fFace(j)%tau0 = (1.+tau0_percentage)*Tdomain%sFault(nf)%fFace(j)%sigma0 * mus
                            else
                                Tdomain%sFault(nf)%fFace(j)%tau0 = -(1.+tau0_percentage)*Tdomain%sFault(nf)%fFace(j)%sigma0 * mus
                            endif
                            Tdomain%sFault(nf)%fFace(j)%mud(:) = mud
                            Tdomain%sFault(nf)%fFace(j)%Ct(:) = Ct
                            Tdomain%sFault(nf)%fFace(j)%b(:) = b
                            Tdomain%sFault(nf)%fFace(j)%w(:) = w
                            Tdomain%sFault(nf)%fFace(j)%Deltau(:,0) = - Deltau_res

                            nv = Tdomain%sFault(nf)%fFace(j)%Face_to_Vertex(0)
                            if (c1*s1*(sigma3-sigma1) > 0.) then
                                Tdomain%sFault(nf)%fVertex(nv)%tau0 =  (1.+tau0_percentage)*Tdomain%sFault(nf)%fVertex(nv)%sigma0 * mus
                            else
                                Tdomain%sFault(nf)%fVertex(nv)%tau0 = - (1.+tau0_percentage)*Tdomain%sFault(nf)%fVertex(nv)%sigma0 * mus
                            endif
                            Tdomain%sFault(nf)%fVertex(nv)%mud = mud
                            Tdomain%sFault(nf)%fVertex(nv)%w = w
                            Tdomain%sFault(nf)%fVertex(nv)%Ct = Ct
                            Tdomain%sFault(nf)%fVertex(nv)%b = b
                            Tdomain%sFault(nf)%fVertex(nv)%Deltau(0) = - Deltau_res
                            nv = Tdomain%sFault(nf)%fFace(j)%Face_to_Vertex(1)
                            if (c1*s1 *(sigma3-sigma1)> 0.) then
                                Tdomain%sFault(nf)%fVertex(nv)%tau0 =  (1.+tau0_percentage)*Tdomain%sFault(nf)%fVertex(nv)%sigma0 * mus
                            else
                                Tdomain%sFault(nf)%fVertex(nv)%tau0 = - (1.+tau0_percentage)*Tdomain%sFault(nf)%fVertex(nv)%sigma0 * mus
                            endif
                            Tdomain%sFault(nf)%fVertex(nv)%mud = mud
                            Tdomain%sFault(nf)%fVertex(nv)%w = w
                            Tdomain%sFault(nf)%fVertex(nv)%Ct = Ct
                            Tdomain%sFault(nf)%fVertex(nv)%b = b
                            Tdomain%sFault(nf)%fVertex(nv)%Deltau(0) = -Deltau_res
                        enddo
                    else
                        do j = 1, 8
                            read (24,*)
                        enddo
                    endif
                enddo
                close (24)
            endif
            nf = nf + 1
        endif
    enddo

    call semname_define_fault_rankl(Tdomain%MPI_Var%my_rank,fnamef)
    open (23,file=fnamef, form="formatted",status="unknown")
    write (23,"(a2)") ">>"
    do n = 0, Tdomain%n_fault-1
        do i = 0, Tdomain%sFault(n)%n_face-1
            nv = Tdomain%sFault(n)%fFace(i)%Face_to_Vertex(0)
            c1 = Tdomain%sFault(n)%Fface(i)%X_Vertex(0)
            write (23,*) c1,Tdomain%sFault(n)%fVertex(nv)%tau0
            do j = 1, ngll-2
                c1 = Tdomain%sFault(n)%Fface(i)%X_Vertex(0) + 0.25 * (Tdomain%sSubDomain(0)%GLLcx(j)+1)
                write (23,*) c1,Tdomain%sFault(n)%fFace(i)%tau0(j)
            enddo
            nv = Tdomain%sFault(n)%fFace(i)%Face_to_Vertex(1)
            c1 = Tdomain%sFault(n)%Fface(i)%X_Vertex(1)
            write (23,*) c1,Tdomain%sFault(n)%fVertex(nv)%tau0
            write (23,"(a2)") ">>"
        enddo
    enddo
    close (23)

    call semname_define_fault_rankn(Tdomain%MPI_Var%my_rank,fnamef)


    open (24,file=fnamef, form="formatted",status="unknown")
    write (24,"(a2)") ">>"
    do n = 0, Tdomain%n_fault-1
        do i = 0, Tdomain%sFault(n)%n_face-1
            nv = Tdomain%sFault(n)%fFace(i)%Face_to_Vertex(0)
            c1 = Tdomain%sFault(n)%Fface(i)%X_Vertex(0)
            write (24,*) c1,Tdomain%sFault(n)%fVertex(nv)%sigma0
            do j = 1, ngll-2
                c1 = Tdomain%sFault(n)%Fface(i)%X_Vertex(0) + 0.25 * (Tdomain%sSubDomain(0)%GLLcx(j)+1)
                write (24,*) c1,Tdomain%sFault(n)%fFace(i)%sigma0(j)
            enddo
            nv = Tdomain%sFault(n)%fFace(i)%Face_to_Vertex(1)
            c1 = Tdomain%sFault(n)%Fface(i)%X_Vertex(1)
            write (24,*) c1,Tdomain%sFault(n)%fVertex(nv)%sigma0
            write (24,"(a2)") ">>"
        enddo
    enddo
    close (24)

    return
end subroutine define_fault_properties
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
