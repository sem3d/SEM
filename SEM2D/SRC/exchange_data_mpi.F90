!>
!!\file exchange_data_MPI.F90
!!\brief this file contains the subroutine that creates the arrays of data
!! that will be exchanged between processors using MPI communications.
!! @author Sebastien terrana
!!\version 1.0
!!\date 21/11/2013
!! This algorithm comes from the following book :
!! cf Algorithm pages 63-64
!<


!>
!! @brief this file contains the subroutine that creates the array of data (forces)
!! that will be sent by the current processor to the processor i_proc 
!! @author Sebastien terrana
!!\version 1.0
!!\date 21/11/2013
!<
subroutine create_send_data (Tdomain,i_proc)
    use constants
  use sdomain
  
  implicit none
  type (domain), intent (INOUT) :: Tdomain
  integer, intent (IN)          :: i_proc

  ! local variables
  integer :: i_send, i_stock, nf, n_face_pointed, ngll, j, nv, nv_aus 

  i_send = Tdomain%Communication_list(i_proc)
  i_stock = 0

  do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
     n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
     ngll = Tdomain%sFace(n_face_pointed)%ngll
     if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
        if(Tdomain%sFace(n_face_pointed)%Type_DG==GALERKIN_CONT) then ! Continuous Galerkin
           Tdomain%sWall(i_proc)%Send_data_2(i_stock:i_stock+ngll-3,0:1) = Tdomain%sFace(n_face_pointed)%Forces (1:ngll-2,0:1)
           i_stock = i_stock + ngll - 2
        else                                              ! Discontinuous Galerkin
           Tdomain%sWall(i_proc)%Send_data_2(i_stock:i_stock+ngll-1,0:4) = Tdomain%sFace(n_face_pointed)%Forces(0:ngll-1,0:4)
           i_stock = i_stock + ngll
        endif
     else
        if(Tdomain%sFace(n_face_pointed)%Type_DG==GALERKIN_CONT) then ! Continuous Galerkin
           do j = 1, ngll-2
              Tdomain%sWall(i_proc)%Send_data_2(i_stock+j-1,0:1) = Tdomain%sFace(n_face_pointed)%Forces(ngll-1-j,0:1)
           enddo
           i_stock = i_stock + ngll - 2
        else                                              ! Discontinuous Galerkin
           do j = 0, ngll-1
              Tdomain%sWall(i_proc)%Send_data_2(i_stock+j-1,0:4) = Tdomain%sFace(n_face_pointed)%Forces(ngll-1-j,0:4)
           enddo
           i_stock = i_stock + ngll
        endif
     endif
  enddo

  ! Treating Vertices (if any)
  do nv = 0, Tdomain%sWall(i_proc)%n_vertices - 1
     nv_aus =  Tdomain%sWall(i_proc)%Vertex_List(nv)
     Tdomain%sWall(i_proc)%Send_data_2 (i_stock,0:1)  = Tdomain%sVertex(nv_aus)%Forces(0:1)
     i_stock = i_stock + 1
  enddo

  ! Treating PML (if any)
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

end subroutine create_send_data


!>
!! @brief this file contains the subroutine that distributes the arrays of data (forces)
!! received by the current processor from the processor i_proc, to the faces at the interfaces
!! between the two processors.
!! @author Sebastien terrana
!!\version 1.0
!!\date 21/11/2013
!<
subroutine assign_recv_data(Tdomain, i_proc)

  use sdomain

  implicit none
  type (domain), intent (INOUT) :: Tdomain
  integer, intent (IN)          :: i_proc

  ! local variables
  integer :: i_send, i_stock, nf, n_face_pointed, ngll, j, nv, nv_aus 

  i_send  = Tdomain%Communication_list(i_proc)
  i_stock = 0

  do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
     n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
     ngll = Tdomain%sFace(n_face_pointed)%ngll
     if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
        if (Tdomain%sFace(n_face_pointed)%Type_DG==GALERKIN_CONT) then ! Continuous Galerkin
           Tdomain%sFace(n_face_pointed)%Forces (1:ngll-2,0:1) =  Tdomain%sFace(n_face_pointed)%Forces ( 1:ngll-2,0:1)  + &
                Tdomain%sWall(i_proc)%Receive_data_2(i_stock:i_stock+ngll-3,0:1)
           i_stock = i_stock + ngll - 2
        else                                               ! Discontinuous Galerkin
           Tdomain%sFace(n_face_pointed)%Forces (0:ngll-1,0:4) =  Tdomain%sFace(n_face_pointed)%Forces (0:ngll-1,0:4)  + &
                Tdomain%sWall(i_proc)%Receive_data_2(i_stock:i_stock+ngll-1,0:4)
           i_stock = i_stock + ngll
        endif
     else
        if (Tdomain%sFace(n_face_pointed)%Type_DG==GALERKIN_CONT) then ! Continuous Galerkin
           do j = 1, ngll-2
              Tdomain%sFace(n_face_pointed)%Forces (ngll-1-j,0:1) =  Tdomain%sFace(n_face_pointed)%Forces (ngll-1-j,0:1)+ &
                   Tdomain%sWall(i_proc)%Receive_data_2(i_stock+j-1,0:1)
           enddo
           i_stock = i_stock + ngll - 2
        else                                               ! Discontinuous Galerkin
           do j = 0, ngll-1
              Tdomain%sFace(n_face_pointed)%Forces (ngll-1-j,0:4) =  Tdomain%sFace(n_face_pointed)%Forces (ngll-1-j,0:4)+ &
                   Tdomain%sWall(i_proc)%Receive_data_2(i_stock+j-1,0:4)
           enddo
           i_stock = i_stock + ngll
        endif
     endif
  enddo

  ! Traiting Vertices (if any)
  do nv = 0, Tdomain%sWall(i_proc)%n_vertices - 1
     nv_aus =  Tdomain%sWall(i_proc)%Vertex_List(nv)
     Tdomain%sVertex(nv_aus)%Forces(0:1) = Tdomain%sVertex(nv_aus)%Forces(0:1) + &
          Tdomain%sWall(i_proc)%Receive_data_2 (i_stock,0:1)
     i_stock = i_stock + 1
  enddo

  ! Traiting PML (if any)
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


end subroutine assign_recv_data
