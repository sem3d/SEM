!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file save_fault_trace.F90
!!\brief Contient la subroutine save_fault_trace.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Assure la sauvegarde des quantites sFault du type domain.
!!
!! \param type (domain), intent (IN) Tdomain
!! \param integer, intent (IN) it
!<


subroutine save_fault_trace (Tdomain, it)

    use sdomain
    use semdatafiles
    use mpi

    implicit none
    type (domain), intent (IN) :: Tdomain
    integer, intent (IN) :: it

    integer   :: i, j, nv, n
    real(fpp) :: dum0, dum1
    character (len=MAX_FILE_SIZE) :: fnamef

    ! local variables
    do n = 0,Tdomain%n_fault-1

#ifdef MKA3D
        call semname_save_fault_trace_rankit("./data/sem/slip_v",Tdomain%MPI_Var%my_rank,it,fnamef)
        open(62,file=fnamef)
        call semname_save_fault_trace_rankit("./data/sem/traction",Tdomain%MPI_Var%my_rank,it,fnamef)
        open(63,file=fnamef)
        call semname_save_fault_trace_rankit("./data/sem/slip_vn",Tdomain%MPI_Var%my_rank,it,fnamef)
        open(64,file=fnamef)
        call semname_save_fault_trace_rankit("./data/sem/tractionn",Tdomain%MPI_Var%my_rank,it,fnamef)
        open(65,file=fnamef)
        call semname_save_fault_trace_rankit("./data/sem/slip",Tdomain%MPI_Var%my_rank,it,fnamef)
        open(66,file=fnamef)

#else
        call semname_save_fault_trace_rankit("slip_v",Tdomain%MPI_Var%my_rank,it,fnamef)
        open(62,file=fnamef)
        call semname_save_fault_trace_rankit("traction",Tdomain%MPI_Var%my_rank,it,fnamef)
        open(63,file=fnamef)
        call semname_save_fault_trace_rankit("slip_vn",Tdomain%MPI_Var%my_rank,it,fnamef)
        open(64,file=fnamef)
        call semname_save_fault_trace_rankit("tractionn",Tdomain%MPI_Var%my_rank,it,fnamef)
        open(65,file=fnamef)
        call semname_save_fault_trace_rankit("slip",Tdomain%MPI_Var%my_rank,it,fnamef)
        open(66,file=fnamef)

#endif

        write (62,"(a2)") ">>"
        write (63,"(a2)") ">>"
        write (64,"(a2)") ">>"
        write (65,"(a2)") ">>"
        write (66,"(a2)") ">>"

        do j = 0, Tdomain%sFault(n)%n_face-1

            nv = Tdomain%sFault(n)%fFace(j)%Face_to_Vertex(0)

            write (62,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(0),&
                Tdomain%sFault(n)%fFace(j)%Z_Vertex(0),&
                Tdomain%sFault(n)%fVertex(nv)%Deltav(0)

            write (63,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(0),&
                Tdomain%sFault(n)%fFace(j)%Z_Vertex(0),&
                Tdomain%sFault(n)%fVertex(nv)%normal(1) *Tdomain%sFault(n)%fVertex(nv)%Traction(0) - &
                Tdomain%sFault(n)%fVertex(nv)%normal(0) * Tdomain%sFault(n)%fVertex(nv)%Traction(1)


            write (64,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(0),&
                Tdomain%sFault(n)%fFace(j)%Z_Vertex(0),&
                Tdomain%sFault(n)%fVertex(nv)%Deltav(1)

            write (65,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(0),&
                Tdomain%sFault(n)%fFace(j)%Z_Vertex(0),&
                Tdomain%sFault(n)%fVertex(nv)%normal(0) * Tdomain%sFault(n)%fVertex(nv)%Traction(0) + &
                Tdomain%sFault(n)%fVertex(nv)%normal(1) * Tdomain%sFault(n)%fVertex(nv)%Traction(1)

            write (66,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(0),&
                Tdomain%sFault(n)%fFace(j)%Z_Vertex(0),&
                Tdomain%sFault(n)%fVertex(nv)%Deltau(0)

            dum0 = 0.5 *  (Tdomain%sFault(n)%fFace(j)%X_Vertex(1)-Tdomain%sFault(n)%fFace(j)%X_Vertex(0))
            dum1 = 0.5 *  (Tdomain%sFault(n)%fFace(j)%Z_Vertex(1)-Tdomain%sFault(n)%fFace(j)%Z_Vertex(0))

            do i = 1, Tdomain%sFault(n)%fFace(j)%ngll-2
                write (62,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(0) + (1. + Tdomain%sSubdomain(0)%GLLcx(i))*dum0, &
                    Tdomain%sFault(n)%fFace(j)%Z_Vertex(0) + (1. + Tdomain%sSubdomain(0)%GLLcx(i))*dum1, &
                    Tdomain%sFault(n)%fFace(j)%Deltav(i,0)
                write (63,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(0) + (1. + Tdomain%sSubdomain(0)%GLLcx(i))*dum0,   &
                    Tdomain%sFault(n)%fFace(j)%Z_Vertex(0) + (1. + Tdomain%sSubdomain(0)%GLLcx(i))*dum1, &
                    Tdomain%sFault(n)%fFace(j)%normal(i,1) * &
                    Tdomain%sFault(n)%fFace(j)%Traction(i,0) -  Tdomain%sFault(n)%fFace(j)%normal(i,0) * &
                    Tdomain%sFault(n)%fFace(j)%Traction(i,1)
                write (64,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(0) + (1. + Tdomain%sSubdomain(0)%GLLcx(i))*dum0,  &
                    Tdomain%sFault(n)%fFace(j)%Z_Vertex(0) + (1. + Tdomain%sSubdomain(0)%GLLcx(i))*dum1, &
                    Tdomain%sFault(n)%fFace(j)%Deltav(i,1)
                write (65,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(0) + (1. + Tdomain%sSubdomain(0)%GLLcx(i))*dum0, &
                    Tdomain%sFault(n)%fFace(j)%Z_Vertex(0) + (1. + Tdomain%sSubdomain(0)%GLLcx(i))*dum1, &
                    Tdomain%sFault(n)%fFace(j)%normal(i,0) * &
                    Tdomain%sFault(n)%fFace(j)%Traction(i,0) +  Tdomain%sFault(n)%fFace(j)%normal(i,1) * &
                    Tdomain%sFault(n)%fFace(j)%Traction(i,1)
                write (66,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(0) + (1. + Tdomain%sSubdomain(0)%GLLcx(i))*dum0,   &
                    Tdomain%sFault(n)%fFace(j)%Z_Vertex(0) + (1. + Tdomain%sSubdomain(0)%GLLcx(i))*dum1, &
                    Tdomain%sFault(n)%fFace(j)%Deltau(i,0)
            enddo

            nv = Tdomain%sFault(n)%fFace(j)%Face_to_Vertex(1)

            write (62,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(1),&
                Tdomain%sFault(n)%fFace(j)%Z_Vertex(1),&
                Tdomain%sFault(n)%fVertex(nv)%Deltav(0)

            write (63,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(1),&
                Tdomain%sFault(n)%fFace(j)%Z_Vertex(1),&
                Tdomain%sFault(n)%fVertex(nv)%normal(1) * Tdomain%sFault(n)%fVertex(nv)%Traction(0) - &
                Tdomain%sFault(n)%fVertex(nv)%normal(0) * Tdomain%sFault(n)%fVertex(nv)%Traction(1)

            write (64,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(1),&
                Tdomain%sFault(n)%fFace(j)%Z_Vertex(1),&
                Tdomain%sFault(n)%fVertex(nv)%Deltav(1)

            write (65,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(1),&
                Tdomain%sFault(n)%fFace(j)%Z_Vertex(1),&
                Tdomain%sFault(n)%fVertex(nv)%normal(0) * Tdomain%sFault(n)%fVertex(nv)%Traction(0) + &
                Tdomain%sFault(n)%fVertex(nv)%normal(1) * Tdomain%sFault(n)%fVertex(nv)%Traction(1)

            write (66,*)Tdomain%sFault(n)%fFace(j)%X_Vertex(1), Tdomain%sFault(n)%fFace(j)%Z_Vertex(1), Tdomain%sFault(n)%fVertex(nv)%Deltau(0)

            write (62,"(a2)") ">>"
            write (63,"(a2)") ">>"
            write (64,"(a2)") ">>"
            write (65,"(a2)") ">>"
            write (66,"(a2)") ">>"
        enddo
    enddo

    call flush(62); close (62); call flush(63); close(63); call flush (64); close (64); call flush (65); close (65)

end subroutine save_fault_trace

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
