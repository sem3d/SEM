!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
subroutine define_Neumann_properties(Tdomain)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer :: ngll,ngll1,ngll2,nf,ne,nv,nf_aus, i_neu, i
    character(len=100) :: surfname
    character(len=12) :: schar

    ! allocations for faces, edges and vertices

    do i=0,Tdomain%n_face-1
      call allocate_face_force(Tdomain%sFace(i))
    enddo
    do i=0,Tdomain%n_edge-1
      call allocate_edge_force(Tdomain%sEdge(i))
    enddo
    do i=0,Tdomain%n_vertex-1
      call allocate_vertex_force(Tdomain%sVertex(i))
    enddo

     do i_neu=lbound(Tdomain%Neumann%NeuSurface,1),ubound(Tdomain%Neumann%NeuSurface,1)
           do nf = 0,Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_faces-1
                ngll1 = Tdomain%Neumann%NeuSurface(i_neu)%Neu_Face(nf)%ngll1
                ngll2 = Tdomain%Neumann%NeuSurface(i_neu)%Neu_Face(nf)%ngll2
                nf_aus = Tdomain%Neumann%NeuSurface(i_neu)%Neu_Face(nf)%Face
                !if(Tdomain%Neumann%Neu_Param%what_bc == "F")then
                !    Tdomain%Neumann%Neu_Face(nf)%mat_index =       &
                !        Tdomain%specel(Tdomain%sFace(nf_aus)%Which_Elem)%mat_index
                !else
                !    Tdomain%Neumann%Neu_Face(nf)%mat_index = Tdomain%Neumann%Neu_Param%mat_index
                !end if
                allocate(Tdomain%Neumann%NeuSurface(i_neu)%Neu_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,0:2))
                allocate(Tdomain%Neumann%NeuSurface(i_neu)%Neu_Face(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2))
                Tdomain%Neumann%NeuSurface(i_neu)%Neu_Face(nf)%BtN = 0d0
            enddo

            do ne = 0,Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_edges-1
                ngll = Tdomain%Neumann%NeuSurface(i_neu)%Neu_Edge(ne)%ngll
                allocate(Tdomain%Neumann%NeuSurface(i_neu)%Neu_Edge(ne)%BtN(1:ngll-2,0:2))
                allocate(Tdomain%Neumann%NeuSurface(i_neu)%Neu_Edge(ne)%Forces(1:ngll-2,0:2))
                Tdomain%Neumann%NeuSurface(i_neu)%Neu_Edge(ne)%BtN = 0d0
            enddo

            do nv = 0,Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_vertices-1
                 Tdomain%Neumann%NeuSurface(i_neu)%Neu_Vertex(nv)%BtN = 0d0
            enddo

            ! elastic properties
            Tdomain%Neumann%Neu_Param%lambda = Tdomain%sSubDomain(Tdomain%Neumann%Neu_Param%mat_index)%DLambda
            Tdomain%Neumann%Neu_Param%Mu = Tdomain%sSubDomain(Tdomain%Neumann%Neu_Param%mat_index)%DMu

            if(Tdomain%Neumann%Neu_Param%wtype == 'S') then
                Tdomain%Neumann%Neu_Param%speed = Tdomain%sSubDomain(Tdomain%Neumann%Neu_Param%mat_index)%Sspeed
            else
                Tdomain%Neumann%Neu_Param%speed = Tdomain%sSubDomain(Tdomain%Neumann%Neu_Param%mat_index)%Pspeed
            endif
       !! print parameters values in the screen
       write(schar,*) i_neu
       surfname = "Neumann"//adjustl(schar(1:len_trim(schar)))
       write(*,*)
       write(*,2004) "SURFACE -> ", i_neu, " :: ", trim(surfname) 
       write(*,2007) "Domaim  : " , Tdomain%Neumann%Neu_Param%mat_index
       write(*,2009) "Lambda  : " , Tdomain%Neumann%Neu_Param%lambda
       write(*,2009) "Mu      : " , Tdomain%Neumann%Neu_Param%Mu
       write(*,2009) "Wave Vp : ", Tdomain%Neumann%Neu_Param%lambda

     enddo   
    write(*,*)
    return
    include 'formats.in'
end subroutine define_Neumann_properties
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine read_neumann_input(Tdomain)
  use sdomain
  
  implicit none
  type(domain), intent(inout)  :: Tdomain
  character(Len=MAX_FILE_SIZE) :: fnamef
  character(len=800)           :: parametric_var
  integer                      :: i, rg
  logical                      :: Boolean(1:20)

  !! Ecriture par un seul proc
  rg = Tdomain%rank
  if (rg.eq.0) then
   WRITE(*,*) 
   WRITE(*,*) "--> READING NEUMANN INPUT"
   WRITE(*,*)  
  endif
  if  (Tdomain%config%neu_Neubc.eq.1) then
       Tdomain%Neumann%Neu_Param%what_bc = 'P'
  elseif (Tdomain%config%neu_Neubc.eq.2) then
       Tdomain%Neumann%Neu_Param%what_bc = 'R'
  elseif (Tdomain%config%neu_Neubc.eq.3) then
       Tdomain%Neumann%Neu_Param%what_bc = 'A'
  elseif (Tdomain%config%neu_Neubc.eq.4) then
       Tdomain%Neumann%Neu_Param%what_bc = 'F'
  elseif (Tdomain%config%neu_Neubc.eq.5) then  
       Tdomain%Neumann%Neu_Param%what_bc = 'S'
  endif

  Tdomain%Neumann%Neu_Param%mat_index = Tdomain%config%neu_mat
  Boolean = .false.
  do i=lbound(Tdomain%config%neu_list,1),ubound(Tdomain%config%neu_list,1)
      Boolean(i) = Tdomain%config%neu_list(i)/=0 
  enddo
  allocate(Tdomain%Neumann%Neu_Param%neu_index(1:Count(Boolean)))
  Tdomain%Neumann%Neu_Param%neu_index = Tdomain%config%neu_list(1:Count(Boolean))
  if (Tdomain%config%neu_type==1) then
      Tdomain%Neumann%Neu_Param%wtype = 'P'
  else
      Tdomain%Neumann%Neu_Param%wtype = 'S'
  end if
  Tdomain%Neumann%Neu_Param%lx = Tdomain%config%neu_L(1)
  Tdomain%Neumann%Neu_Param%ly = Tdomain%config%neu_L(2)
  Tdomain%Neumann%Neu_Param%lz = Tdomain%config%neu_L(3)
  Tdomain%Neumann%Neu_Param%xs = Tdomain%config%neu_C(1)
  Tdomain%Neumann%Neu_Param%ys = Tdomain%config%neu_C(2)
  Tdomain%Neumann%Neu_Param%zs = Tdomain%config%neu_C(3)
  Tdomain%Neumann%Neu_Param%f0 = Tdomain%config%neu_f0

  !! Add by Mtaro
  if (Tdomain%Neumann%Neu_Param%what_bc=="A") then
      Tdomain%Neumann%Neu_Param%neu_dim   = Tdomain%config%neu_dim
      Tdomain%Neumann%Neu_Param%neu_varia(1:LEN_TRIM(fromcstr(Tdomain%config%neu_varia))) = &
                                            TRIM(fromcstr(Tdomain%config%neu_varia))
      Tdomain%Neumann%Neu_Param%neu_source(1:LEN_TRIM(fromcstr(Tdomain%config%neu_source)))= &
                                             TRIM(fromcstr(Tdomain%config%neu_source))
      Tdomain%Neumann%Neu_Param%neu_funcx(1:LEN_TRIM(fromcstr(Tdomain%config%neu_funcx))) = &
                                            TRIM(fromcstr(Tdomain%config%neu_funcx))
      if (Tdomain%Neumann%Neu_Param%neu_dim.gt.1) then
          Tdomain%Neumann%Neu_Param%neu_funcy(1:LEN_TRIM(fromcstr(Tdomain%config%neu_funcy))) = &
                                                TRIM(fromcstr(Tdomain%config%neu_funcy))
      endif
      if (Tdomain%Neumann%Neu_Param%neu_dim.gt.2) then
          Tdomain%Neumann%Neu_Param%neu_funcz(1:LEN_TRIM(fromcstr(Tdomain%config%neu_funcz))) = &
                                                TRIM(fromcstr(Tdomain%config%neu_funcz))
      endif
      if (Tdomain%Neumann%Neu_Param%neu_source.eq."M") then
           Tdomain%Neumann%Neu_Param%neu_funcxy(1:LEN_TRIM(fromcstr(Tdomain%config%neu_funcxy)))= &
                                                  TRIM(fromcstr(Tdomain%config%neu_funcxy))
           if (Tdomain%Neumann%Neu_Param%neu_dim.eq.3) then
              Tdomain%Neumann%Neu_Param%neu_funcxz(1:LEN_TRIM(fromcstr(Tdomain%config%neu_funcxz)))= &
                                                     TRIM(fromcstr(Tdomain%config%neu_funcxz))
              Tdomain%Neumann%Neu_Param%neu_funcyz(1:LEN_TRIM(fromcstr(Tdomain%config%neu_funcyz)))= &
                                                    TRIM(fromcstr(Tdomain%config%neu_funcyz))
           endif
      endif
     !!
     !!
      Tdomain%Neumann%Neu_Param%neu_paramvar=Tdomain%config%neu_paramvar
      if (Tdomain%Neumann%Neu_Param%neu_paramvar==1) then
         Tdomain%Neumann%Neu_Param%neu_nparamvar = Tdomain%config%neu_nparamvar
         allocate(Tdomain%Neumann%Neu_Param%neu_paravalue(1:Tdomain%config%neu_nparamvar))
         do i=1,Tdomain%config%neu_nparamvar
           Tdomain%Neumann%Neu_Param%neu_paravalue(i) = Tdomain%config%neu_Paravalue(i)
         enddo
         parametric_var(1:LEN_TRIM(fromcstr(Tdomain%config%neu_paramname)))=TRIM(fromcstr(Tdomain%config%neu_paramname))
         call Split(parametric_var, Tdomain%config%neu_nparamvar, &
                      Tdomain%Neumann%Neu_Param%neu_paramname)
      endif
      !!
      !!
      if (rg.eq.0) then
      write(*,*)
      write(*,*) " Neuman BC Source analytical defined "
      write(*,*) "-------------------------------------"
      write(*,1006) trim("Va = "//adjustr(Tdomain%Neumann%Neu_Param%neu_varia))
      write(*,1005) trim("F1 = "//adjustr(Tdomain%Neumann%Neu_Param%neu_funcx))
      if (Tdomain%Neumann%Neu_Param%neu_dim.gt.1) then
         write(*,1005) trim("F2 = "//adjustr(Tdomain%Neumann%Neu_Param%neu_funcy))
      endif
      if (Tdomain%Neumann%Neu_Param%neu_dim.gt.2) then
         write(*,1005) trim("F3 = "//adjustr(Tdomain%Neumann%Neu_Param%neu_funcz))
      endif
      if (Tdomain%Neumann%Neu_Param%neu_source.eq."M") then
        write(*,1005) trim("F4 = "//adjustr(Tdomain%Neumann%Neu_Param%neu_funcxy))
        if (Tdomain%Neumann%Neu_Param%neu_dim.gt.2) then
           write(*,1005) trim("F5 = "//adjustr(Tdomain%Neumann%Neu_Param%neu_funcxz))
           write(*,1005) trim("F6 = "//adjustr(Tdomain%Neumann%Neu_Param%neu_funcyz))
        endif
     endif
     if (Tdomain%Neumann%Neu_Param%neu_paramvar==1) then
       write(*,1005) trim("Param ="//adjustr(parametric_var))
       write(*,2014) "ParamVal =",Tdomain%Neumann%Neu_Param%neu_paravalue
     endif
     endif
     write(*,*)
 endif
  
  include 'formats.in'
end subroutine read_neumann_input

subroutine Split(String, nstring, SubStrings)
      implicit none
      
      character(len=*), intent(in)                         :: String
      integer, intent(in)                                  :: nstring
      character, parameter                                 :: Delimiter= " "
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
