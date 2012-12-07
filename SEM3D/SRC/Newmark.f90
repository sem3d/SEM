subroutine Newmark(Tdomain,rg,ntime)
 ! Predictor-MultiCorrector Newmark Velocity Scheme within a 
 ! Time staggered Stress-Velocity formulation inside PML
    use sdomain
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: rg,ntime

    logical, dimension(:), allocatable :: L_Face, L_Edge, L_Vertex
    integer ::  n, mat, nelem, ngll1,ngll2,ngll3, ngllx,nglly,ngllz,    &
                nf,ne,nv, nf_aus,ne_aus,nv_aus,ipoint,ndt2, ngll, code, &
                i,j,k, x,y,z,ngllPML,shift, I_give_to, I_take_from,     &
                n_rings, ntimetrace,flag,nummax,ngll_F,ngllPML_F
    integer, parameter :: etiquette = 100
    integer, dimension(mpi_status_size) :: statut
    real :: alpha, bega, gam1, Dt, max, maximum,t1,t2,time_dep
    real, dimension(0:2) :: V_free_vertex
    real, dimension(:,:), allocatable :: VE_free
    real, dimension(:,:,:), allocatable :: VF_free
    logical :: sortie



if(.not. Tdomain%TimeD%velocity_scheme)   &
         stop "Newmark scheme implemented only in velocity form."

!- Prediction Phase   
call Newmark_Predictor(Tdomain,rg)

!- Solution phase
call internal_forces(Tdomain,rg)

! External Forces
if(Tdomain%logicD%any_source)then
    call external_forces(Tdomain,Tdomain%TimeD%rtime,rg)
end if

! Communication of Forces within a single process
call inside_proc_forces(Tdomain)

! MPI communications
if(Tdomain%n_proc > 1)then
  ! from external faces, edges and vertices to Communication global arrays
    do n = 0,Tdomain%n_proc-1
        call Comm_Forces_Complete(n,Tdomain)
        call Comm_Forces_PML_Complete(n,Tdomain)
    end do

  ! now we can exchange force values with proc n
    n = Tdomain%n_proc
    do shift = 1,n-1
        call shift_to_parameters(rg,n,shift,I_give_to,I_take_from,n_rings)
    ! solid forces exchange
        call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngll,        &
               Tdomain%sComm(I_take_from)%ngll,I_give_to,I_take_from,3,   &
               Tdomain%sComm(I_give_to)%GiveForces,Tdomain%sComm(I_take_from)%TakeForces)
        call MPI_BARRIER(MPI_COMM_WORLD,code)
    ! fluid forces exchange
        call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngll_F,        &
               Tdomain%sComm(I_take_from)%ngll_F,I_give_to,I_take_from,1,   &
               Tdomain%sComm(I_give_to)%GiveForcesFl,Tdomain%sComm(I_take_from)%TakeForcesFl)
        call MPI_BARRIER(MPI_COMM_WORLD,code)
    ! solid PML forces exchange
        call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngllPML,        &
               Tdomain%sComm(I_take_from)%ngllPML,I_give_to,I_take_from,9,   &
               Tdomain%sComm(I_give_to)%GiveForcesPML,Tdomain%sComm(I_take_from)%TakeForcesPML)
        call MPI_BARRIER(MPI_COMM_WORLD,code)
    ! fluid PML forces exchange
        call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngllPML_F,        &
               Tdomain%sComm(I_take_from)%ngllPML_F,I_give_to,I_take_from,3,   &
               Tdomain%sComm(I_give_to)%GiveForcesPMLFl,Tdomain%sComm(I_take_from)%TakeForcesPMLFl)
        call MPI_BARRIER(MPI_COMM_WORLD,code)

    end do  ! do shift

  ! now: assemblage on external faces, edges and vertices
    do n = 0,Tdomain%n_proc-1
        ngll = 0
        ngll_F = 0
        ngllPML = 0
        ngllPML_F = 0
        call Comm_Forces_Face(Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
        call Comm_Forces_Edge(Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
        call Comm_Forces_Vertex(Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
    enddo

endif  ! if nproc > 1


! Neumann B.C.: associated forces
if(Tdomain%logicD%neumann_local_present)then
   do nf = 0,Tdomain%Neumann%Neu_n_faces-1
      ngll1 = Tdomain%Neumann%Neu_Face(nf)%ngll1
      ngll2 = Tdomain%Neumann%Neu_Face(nf)%ngll2
      mat = Tdomain%Neumann%Neu_Face(nf)%mat_index
      call compute_Neu_forces_on_face(Tdomain%Neumann%Neu_Face(nf),     &
             Tdomain%Neumann%Neu_Param,Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime)
   enddo
   do ne = 0,Tdomain%Neumann%Neu_n_edges-1
      ngll = Tdomain%Neumann%Neu_Edge(ne)%ngll
      ne_aus = Tdomain%Neumann%Neu_Edge(ne)%Edge
      mat = Tdomain%sEdge(ne_aus)%mat_index
      call compute_Neu_forces_on_edge(Tdomain%Neumann%Neu_Edge(ne),     &
             Tdomain%Neumann%Neu_Param,Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime)
   enddo
   do nv = 0,Tdomain%Neumann%Neu_n_vertices-1   
      nv_aus = Tdomain%Neumann%Neu_Vertex(nv)%Vertex
      mat = Tdomain%sVertex(nv_aus)%mat_index
      n = merge(0,1,nv == 4)
      call compute_Neu_forces_on_vertex(Tdomain%Neumann%Neu_Vertex(nv),n,  &
             Tdomain%Neumann%Neu_Param,Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime)
   enddo
  ! addition of Neumann forces
   do nf = 0, Tdomain%Neumann%Neu_n_faces-1
      nf_aus = Tdomain%Neumann%Neu_Face(nf)%Face
      if(.not.Tdomain%sface(nf_aus)%PML)then
        Tdomain%sFace(nf_aus)%Forces(:,:,0:2) = Tdomain%sFace(nf_aus)%Forces(:,:,0:2) - & 
                                             Tdomain%Neumann%Neu_Face(nf)%Forces(:,:,0:2)
      else
        Tdomain%sFace(nf_aus)%Forces3(:,:,0:2) = Tdomain%sFace(nf_aus)%Forces3(:,:,0:2) - & 
                                             Tdomain%Neumann%Neu_Face(nf)%Forces(:,:,0:2)
     endif
   enddo
   do ne = 0, Tdomain%Neumann%Neu_n_edges-1
      ne_aus = Tdomain%Neumann%Neu_Edge(ne)%Edge
      if(.not.Tdomain%sedge(ne_aus)%PML)then
        Tdomain%sEdge(ne_aus)%Forces(:,0:2) = Tdomain%sEdge(ne_aus)%Forces(:,0:2) - & 
                                            Tdomain%Neumann%Neu_Edge(ne)%Forces(:,0:2)
      else
        Tdomain%sEdge(ne_aus)%Forces3(:,0:2) = Tdomain%sEdge(ne_aus)%Forces3(:,0:2) - & 
                                            Tdomain%Neumann%Neu_Edge(ne)%Forces(:,0:2)
      endif
   enddo
   do nv = 0, Tdomain%Neumann%Neu_n_vertices-1
      nv_aus = Tdomain%Neumann%Neu_Vertex(nv)%Vertex
      if(.not.Tdomain%svertex(nv_aus)%PML)then
        Tdomain%sVertex(nv_aus)%Forces(0:2) = Tdomain%sVertex(nv_aus)%Forces(0:2) -  &
                                              Tdomain%Neumann%Neu_Vertex(nv)%Forces(0:2)
      else
        Tdomain%sVertex(nv_aus)%Forces3(0:2) = Tdomain%sVertex(nv_aus)%Forces3(0:2) - &
                                               Tdomain%Neumann%Neu_Vertex(nv)%Forces(0:2)
      endif
   enddo
endif


! Correction phase
allocate (L_Face(0:Tdomain%n_face-1))
L_Face = .true.
allocate (L_Edge(0:Tdomain%n_edge-1))
L_Edge = .true.
allocate (L_Vertex(0:Tdomain%n_vertex-1))
L_Vertex = .true.
do n = 0,Tdomain%n_elem-1
    mat = Tdomain%specel(n)%mat_index
    if (.not. Tdomain%specel(n)%PML) then
        call Correction_Elem_Veloc (Tdomain%specel(n),bega, gam1,Tdomain%sSubDomain(mat)%Dt)
    else 
        if (Tdomain%specel(n)%FPML) then
             call Correction_Elem_FPML_Veloc (Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt, Tdomain%sSubdomain(mat)%freq)
        else
            call Correction_Elem_PML_Veloc (Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt)
        endif
    endif
    do i = 0,5
        nf = Tdomain%specel(n)%Near_Faces(i)
        if (L_Face(nf)) then
            L_Face(nf) = .false.
            if (.not.Tdomain%sface(nf)%PML) then
                call Correction_Face_Veloc (Tdomain%sface(nf), bega, gam1, dt)
            else
              if (Tdomain%sface(nf)%FPML) then
                   call Correction_Face_FPML_Veloc (Tdomain%sface(nf), dt, Tdomain%sSubdomain(mat)%freq)
              else
                   call Correction_Face_PML_Veloc (Tdomain%sface(nf), dt)
              endif
            endif
        endif
    enddo
    do i = 0,11
        ne = Tdomain%specel(n)%Near_Edges(i)
        if (L_Edge(ne)) then
            L_Edge(ne) = .false.
            if (.not.Tdomain%sedge(ne)%PML) then
                call Correction_Edge_Veloc (Tdomain%sedge(ne), bega, gam1, dt,ne)
            else
              if (Tdomain%sedge(ne)%FPML) then
                   call Correction_Edge_FPML_Veloc (Tdomain%sedge(ne), dt, Tdomain%sSubdomain(mat)%freq)
             else
                   call Correction_Edge_PML_Veloc (Tdomain%sedge(ne), dt)
             endif
            endif
        endif
    enddo
    do i = 0,7
        nv = Tdomain%specel(n)%Near_Vertices(i)
        if(L_Vertex(nv))then
            L_Vertex(nv) = .false.
            if(.not. Tdomain%svertex(nv)%PML)then
                call Correction_Vertex_Veloc (Tdomain%svertex(nv), bega, gam1, dt)
            else
              if(Tdomain%svertex(nv)%FPML)then
                  call Correction_Vertex_FPML_Veloc(Tdomain%svertex(nv),dt,Tdomain%sSubdomain(mat)%freq)
              else
                  call Correction_Vertex_PML_Veloc(Tdomain%svertex(nv),dt)
             endif
            endif
        endif
    enddo
enddo
deallocate (L_Face,L_Edge,L_Vertex)


! Save Trace

if (Tdomain%logicD%save_trace) then

   do n = 0, Tdomain%n_receivers-1

      ndt2 = Tdomain%sReceiver(n)%ndt

      if (rg == Tdomain%sReceiver(n)%proc) then

         if ( mod (ntime,Tdomain%TimeD%ntrace) == 0 ) then
          if ( Tdomain%sReceiver(n)%flag == 1 ) allocate (Tdomain%sReceiver(n)%StoreTrace (0:Tdomain%TimeD%ntrace-1, 0:2) )
          if ( Tdomain%sReceiver(n)%flag == 2 ) allocate (Tdomain%sReceiver(n)%StoreTrace (0:(Tdomain%TimeD%ntrace-1)/ndt2, 0:2) )
          Tdomain%sReceiver(n)%StoreTrace = 0.
        endif

        i = Tdomain%sReceiver(n)%elem
        ngll1 = Tdomain%specel(i)%ngllx
        ngll2 = Tdomain%specel(i)%nglly
        ngll3 = Tdomain%specel(i)%ngllz

        if ( Tdomain%sReceiver(n)%flag == 1 .or. ( (Tdomain%sReceiver(n)%flag == 2) .and. (mod (ntime+1,ndt2)==0) ) ) then    

         if ( Tdomain%sReceiver(n)%flag == 1 ) ntimetrace = mod (ntime,Tdomain%TimeD%ntrace)
         if ( Tdomain%sReceiver(n)%flag == 2 ) ntimetrace = mod ((ntime+1)/ndt2-1,Tdomain%TimeD%ntrace/ndt2)

         do x = 0,ngll1-1
          do y = 0,ngll2-1
           do z = 0,ngll3-1
              if (x==0) then
                 if (y==0) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(0)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(4)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(6)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(z,:)
                    endif
                 else if (y==ngll2-1) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(3)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(7)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(10)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(z,:)
                    endif
                 else if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(3)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(y,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(11)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(y,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(4)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(y,z,:)
                 endif
              else if (x==ngll1-1) then
                 if (y==0) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(1)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(5)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(4)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(z,:)
                    endif
                 else if (y==ngll2-1) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(2)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(6)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(7)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(z,:)
                    endif
                 else if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(1)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(y,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(8)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(y,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(2)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(y,z,:)
                 endif
              else if (y==0) then
                 if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(0)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(x,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(5)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(x,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(1)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(x,z,:)
                 endif
              else if (y==ngll2-1) then
                 if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(2)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(x,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(9)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(x,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(3)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(x,z,:)
                 endif
              else if (z==0) then
                 nf = Tdomain%specel(i)%Near_Faces(0)
                 Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(x,y,:)
              else if (z==ngll3-1) then
                 nf = Tdomain%specel(i)%Near_Faces(5)
                 Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(x,y,:)
              else
                 Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%specel(i)%Veloc(x,y,z,:)
              endif
           enddo
          enddo
         enddo

         do x = 0,ngll1-1
          do y = 0,ngll2-1
           do z = 0,ngll3-1
              Tdomain%sReceiver(n)%StoreTrace(ntimetrace,:) = Tdomain%sReceiver(n)%StoreTrace(ntimetrace,:) + &
                                Tdomain%sReceiver(n)%coeff(x,y,z,:) * Tdomain%sReceiver(n)%pol(x,y,z)
           enddo
          enddo
         enddo

       endif

      endif
   enddo
endif


if(rg == 0) print *,'Num iteration - Time',ntime,Tdomain%TimeD%rtime

return

end subroutine Newmark
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine Newmark_Predictor(Tdomain,rg)

    use sdomain
    implicit none
    
    type(domain), intent(inout)   :: Tdomain
    integer, intent(in)  :: rg
    real  :: alpha,bega,gam1,dt
    integer  :: i,n,mat,nf,ne,nv
    logical, dimension(:), allocatable  :: L_Face,L_Edge,L_Vertex

alpha = Tdomain%TimeD%alpha
bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
gam1 = 1. / Tdomain%TimeD%gamma


do n = 0,Tdomain%n_elem-1 
    mat = Tdomain%specel(n)%mat_index 
    dt = Tdomain%sSubDomain(mat)%dt
   ! SOLID PART
    if(Tdomain%specel(n)%solid)then
        if(.not. Tdomain%specel(n)%PML)then  ! physical part
            call Prediction_Elem_Veloc(Tdomain%specel(n),alpha,bega,gam1,Dt)
        else    ! PML part
            call get_PMLprediction_v2el(Tdomain,n,bega,dt,rg) 
            call get_PMLprediction_e2el(Tdomain,n,bega,dt,rg)
            call get_PMLprediction_f2el(Tdomain,n,bega,dt,rg)
            if(Tdomain%specel(n)%FPML)then
                call Prediction_Elem_FPML_Veloc(Tdomain%specel(n),bega,dt,               &
                      Tdomain%sSubDomain(mat)%hTPrimex,Tdomain%sSubDomain(mat)%hPrimey,  &
                      Tdomain%sSubDomain(mat)%hprimez,rg,n,Tdomain%sSubDomain(mat)%freq)
            else
                call Prediction_Elem_PML_Veloc(Tdomain%specel(n),bega,dt,                &
                      Tdomain%sSubDomain(mat)%hTPrimex,Tdomain%sSubDomain(mat)%hPrimey,  &
                      Tdomain%sSubDomain(mat)%hprimez,rg,n)
            endif
        endif
   ! FLUID PART
    else   
        if(.not. Tdomain%specel(n)%PML)then  ! physical part
            call Prediction_Elem_VelPhi(Tdomain%specel(n),alpha,bega,gam1,Dt)
        else    ! PML part
            call get_PMLprediction_v2el_fl(Tdomain,n,bega,dt,rg) 
            call get_PMLprediction_e2el_fl(Tdomain,n,bega,dt,rg)
            call get_PMLprediction_f2el_fl(Tdomain,n,bega,dt,rg)
            call Prediction_Elem_PML_VelPhi(Tdomain%specel(n),bega,dt,             &
                 Tdomain%sSubDomain(mat)%hTPrimex,Tdomain%sSubDomain(mat)%hPrimey, &
                 Tdomain%sSubDomain(mat)%hprimez)
        endif
    end if
enddo

allocate(L_Face(0:Tdomain%n_face-1))
L_Face = .true.
allocate(L_Edge(0:Tdomain%n_edge-1))
L_Edge = .true.
allocate(L_Vertex(0:Tdomain%n_vertex-1))
L_Vertex = .true.
do n = 0,Tdomain%n_elem-1
    mat = Tdomain%specel(n)%mat_index
    dt = Tdomain%sSubDomain(mat)%dt
    do i = 0,5
        nf = Tdomain%specel(n)%Near_Faces(i)
        if(L_Face(nf))then
            L_Face(nf) = .false.
            if(.not.Tdomain%sface(nf)%PML)then
                if(Tdomain%sFace(nf)%solid)then
                    call Prediction_Face_Veloc(Tdomain%sface(nf),alpha,bega,gam1,dt)
                else
                    call Prediction_Face_VelPhi(Tdomain%sface(nf),alpha,bega,gam1,dt)
                end if
            end if
        endif
    enddo
    do i = 0,11
        ne = Tdomain%specel(n)%Near_Edges(i)
        if(L_Edge(ne))then
            L_Edge(ne) = .false.
            if(.not.Tdomain%sedge(ne)%PML)then
                if(Tdomain%sEdge(ne)%solid)then
                    call Prediction_Edge_Veloc(Tdomain%sedge(ne),alpha,bega,gam1,dt)
                else
                    call Prediction_Edge_VelPhi(Tdomain%sedge(ne),alpha,bega,gam1,dt)
                end if
            end if
        endif
    enddo
    do i = 0,7
        nv = Tdomain%specel(n)%Near_Vertices(i)
        if(L_Vertex(nv))then
            L_Vertex(nv) = .false.
            if(.not.Tdomain%svertex(nv)%PML)then
                if(Tdomain%sVertex(nv)%solid)then
                    call Prediction_Vertex_Veloc(Tdomain%svertex(nv),alpha,bega,gam1,dt)
                else
                    call Prediction_Vertex_VelPhi(Tdomain%svertex(nv),alpha,bega,gam1,dt)
                end if
            end if
        endif
    enddo
enddo
deallocate(L_Face,L_Edge,L_Vertex)

return

end subroutine Newmark_Predictor 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine internal_forces(Tdomain,rank)
 ! volume forces - depending on rheology
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: rank
    integer  :: n,mat


do n = 0,Tdomain%n_elem-1
    mat = Tdomain%specel(n)%mat_index
    if(Tdomain%specel(n)%solid)then  ! solid part
        if(.not.Tdomain%specel(n)%PML)then
            call get_Displ_Face2Elem(Tdomain,n)
            call get_Displ_Edge2Elem(Tdomain,n)
            call get_Displ_Vertex2Elem(Tdomain,n)
            call compute_InternalForces_Elem(Tdomain%specel(n),                    &
                 Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimex, &
                 Tdomain%sSubDomain(mat)%hprimey,Tdomain%sSubDomain(mat)%hTprimey, &
                 Tdomain%sSubDomain(mat)%hprimez,Tdomain%sSubDomain(mat)%hTprimez)
        else
            call compute_InternalForces_PML_Elem(Tdomain%specel(n),                &
                 Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimey, &
                 Tdomain%sSubDomain(mat)%hTprimez)
        end if
    else   ! fluid part
        if(.not.Tdomain%specel(n)%PML)then
            call get_Phi_Face2Elem(Tdomain,n)
            call get_Phi_Edge2Elem(Tdomain,n)
            call get_Phi_Vertex2Elem(Tdomain,n)
            call compute_InternalForces_Elem_Fluid(Tdomain%specel(n),              &
                 Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimex, &
                 Tdomain%sSubDomain(mat)%hprimey,Tdomain%sSubDomain(mat)%hTprimey, &
                 Tdomain%sSubDomain(mat)%hprimez,Tdomain%sSubDomain(mat)%hTprimez)
        else
            call compute_InternalForces_PML_Elem_Fl(Tdomain%specel(n),             &
                 Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimey, &
                 Tdomain%sSubDomain(mat)%hTprimez)
        end if

    end if
end do

return
end subroutine internal_forces
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine external_forces(Tdomain,timer,rank)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)  :: rank
    real, intent(in)  :: timer
    integer  :: ns,nel,ngllx,nglly,ngllz,i_dir
    real :: t

do ns = 0, Tdomain%n_source-1
   if(rank == Tdomain%sSource(ns)%proc)then
       nel = Tdomain%Ssource(ns)%elem
       ngllx = Tdomain%specel(nel)%ngllx
       nglly = Tdomain%specel(nel)%nglly
       ngllz = Tdomain%specel(nel)%ngllz
     ! time : t_(n+1/2) for solid ; t_n for fluid 
       t = merge(timer+Tdomain%TimeD%dtmin/2d0,timer,Tdomain%specel(nel)%solid)
     !
       if(Tdomain%sSource(ns)%i_type_source == 1)then  ! collocated force in solid
           i_dir = Tdomain%Ssource(ns)%i_dir
           Tdomain%specel(nel)%Forces(:,:,:,i_dir) = Tdomain%specel(nel)%Forces(:,:,:,i_dir)+ &
                 CompSource(Tdomain%sSource(ns),t,i_dir)*Tdomain%sSource(ns)%ExtForce(:,:,:)
       else if(Tdomain%sSource(ns)%i_type_source == 2)then  ! moment tensor source
         do i_dir = 0,2
           Tdomain%specel(nel)%Forces(:,:,:,i_dir) = Tdomain%specel(nel)%Forces(:,:,:,i_dir)+ &
                 CompSource(Tdomain%sSource(ns),t,i_dir)*Tdomain%sSource(ns)%MomForce(:,:,:,i_dir)
         end do
       else if(Tdomain%sSource(ns)%i_type_source == 3)then    ! pressure pulse in fluid
           Tdomain%specel(nel)%ForcesFl(:,:,:) = Tdomain%specel(nel)%ForcesFl(:,:,:)+    &
                 CompSource_Fl(Tdomain%sSource(ns),t)*Tdomain%sSource(ns)%ExtForce(:,:,:)
       end if
   endif
enddo

return
end subroutine external_forces
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine inside_proc_forces(Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer  :: n,nf,ne,nv

! init.
do nf = 0,Tdomain%n_face-1
    if(Tdomain%sFace(nf)%solid)then
        Tdomain%sFace(nf)%Forces = 0
        if(Tdomain%sFace(nf)%PML)then
            Tdomain%sFace(nf)%Forces1 = 0
            Tdomain%sFace(nf)%Forces2 = 0
            Tdomain%sFace(nf)%Forces3 = 0
        endif
    else
        Tdomain%sFace(nf)%ForcesFl = 0
        if(Tdomain%sFace(nf)%PML)then
            Tdomain%sFace(nf)%ForcesFl1 = 0
            Tdomain%sFace(nf)%ForcesFl2 = 0
            Tdomain%sFace(nf)%ForcesFl3 = 0
        endif
    end if
enddo
do ne = 0,Tdomain%n_edge-1
    if(Tdomain%sEdge(ne)%solid)then
        Tdomain%sEdge(ne)%Forces = 0
        if(Tdomain%sEdge(ne)%PML)then
            Tdomain%sEdge(ne)%Forces1 = 0
            Tdomain%sEdge(ne)%Forces2 = 0
            Tdomain%sEdge(ne)%Forces3 = 0
        endif
    else
        Tdomain%sEdge(ne)%ForcesFl = 0
        if(Tdomain%sEdge(ne)%PML)then
            Tdomain%sEdge(ne)%ForcesFl1 = 0
            Tdomain%sEdge(ne)%ForcesFl2 = 0
            Tdomain%sEdge(ne)%ForcesFl3 = 0
        end if
    end if
enddo
do nv = 0,Tdomain%n_vertex-1
    if(Tdomain%sVertex(nv)%solid)then
        Tdomain%sVertex(nv)%Forces = 0
        if(Tdomain%sVertex(nv)%PML)then
            Tdomain%sVertex(nv)%Forces1 = 0
            Tdomain%sVertex(nv)%Forces2 = 0
            Tdomain%sVertex(nv)%Forces3 = 0
        endif
    else
        Tdomain%sVertex(nv)%ForcesFl = 0
        if(Tdomain%sVertex(nv)%PML)then
            Tdomain%sVertex(nv)%ForcesFl1 = 0
            Tdomain%sVertex(nv)%ForcesFl2 = 0
            Tdomain%sVertex(nv)%ForcesFl3 = 0
        endif

    end if
enddo

do n = 0,Tdomain%n_elem-1
    if(Tdomain%specel(n)%solid)then
        call get_Forces_Elem2Face(Tdomain,n)
        call get_Forces_Elem2Edge(Tdomain,n)
        call get_Forces_Elem2Vertex(Tdomain,n)
    else
        call get_ForcesFl_Elem2Face(Tdomain,n)
        call get_ForcesFl_Elem2Edge(Tdomain,n)
        call get_ForcesFl_Elem2Vertex(Tdomain,n)
    end if
enddo

return
end subroutine inside_proc_forces
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine Comm_Forces_Complete(n,Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: n
    integer  :: ngll,ngll_F,i,j,k,nf,ne,nv

    ngll = 0 ; ngll_F = 0
  ! faces
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        if(Tdomain%sFace(nf)%solid)then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sFace(nf)%Forces(k,j,0:2)
                    ngll = ngll + 1
                enddo
            enddo
        else
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GiveForcesFl(ngll_F) = Tdomain%sFace(nf)%ForcesFl(k,j)
                    ngll_F = ngll_F + 1
                enddo
            enddo
        end if
    enddo
  ! edges
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        if(Tdomain%sEdge(ne)%solid)then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2)
                ngll = ngll + 1
            enddo
        else
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%GiveForcesFl(ngll_F) = Tdomain%sEdge(ne)%ForcesFl(j)
                ngll_F = ngll_F + 1
            enddo
        end if
    enddo
  ! vertices
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        if(Tdomain%sVertex(nv)%solid)then
            Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sVertex(nv)%Forces(0:2)
            ngll = ngll + 1
        else
            Tdomain%sComm(n)%GiveForcesFl(ngll_F) = Tdomain%sVertex(nv)%ForcesFl
            ngll_F = ngll_F + 1
        end if
    enddo

return
end subroutine Comm_Forces_Complete
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine Comm_Forces_PML_Complete(n,Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: n
    integer  :: ngllPML,ngllPML_F,i,j,k,nf,ne,nv

    ngllPML = 0 ; ngllPML_F = 0

  ! faces
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        if(Tdomain%sFace(nf)%PML)then
          if(Tdomain%sFace(nf)%solid)then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sFace(nf)%Forces1(k,j,0:2)
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sFace(nf)%Forces2(k,j,0:2)
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sFace(nf)%Forces3(k,j,0:2)
                    ngllPML = ngllPML + 1
                enddo
            enddo
          else
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,1) = Tdomain%sFace(nf)%ForcesFl1(k,j)
                    Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,2) = Tdomain%sFace(nf)%ForcesFl2(k,j)
                    Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,3) = Tdomain%sFace(nf)%ForcesFl3(k,j)
                    ngllPML_F = ngllPML_F + 1
                enddo
            enddo
          end if
        endif
    enddo
  ! edges
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        if(Tdomain%sEdge(ne)%PML)then
          if(Tdomain%sEdge(ne)%solid)then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2)
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2)
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2)
                ngllPML = ngllPML + 1
            enddo
          else
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,1) = Tdomain%sEdge(ne)%ForcesFl1(j)
                Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,2) = Tdomain%sEdge(ne)%ForcesFl2(j)
                Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,3) = Tdomain%sEdge(ne)%ForcesFl3(j)
                ngllPML_F = ngllPML_F + 1
            enddo
          end if
        endif
    enddo
  ! vertices
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        if(Tdomain%sVertex(nv)%PML)then
          if(Tdomain%sVertex(nv)%solid)then
            Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sVertex(nv)%Forces1(0:2)
            Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sVertex(nv)%Forces2(0:2)
            Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sVertex(nv)%Forces3(0:2)
            ngllPML = ngllPML + 1
          else
            Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,1) = Tdomain%sVertex(nv)%ForcesFl1
            Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,2) = Tdomain%sVertex(nv)%ForcesFl2
            Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,3) = Tdomain%sVertex(nv)%ForcesFl3
            ngllPML_F = ngllPML_F + 1
          end if
        endif
    enddo

return
end subroutine Comm_Forces_PML_Complete
