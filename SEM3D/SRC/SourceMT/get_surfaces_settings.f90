Module surfSetting

use sdomain
use sem_hdf5
use mpi
use sem_c_bindings
use semdatafiles, only : MAX_FILE_SIZE          
use constants 

contains

        subroutine read_neu_surface(Tdomain)

         implicit none
         !
         type(domain), intent(inout)          :: Tdomain
         !
         integer(HID_T)                       :: fid, proc_id, surf_id 
         integer                              :: hdferr, ierr, i_neu, i_surf, i, index_neu
         integer, allocatable, dimension(:,:) :: itemp2, itemp2b, itemp2c
         integer, allocatable, dimension(:,:) :: itemp
         character(len=12) ::  schar
         character(len=100)                   :: surfacename

           do i_neu=0,size(Tdomain%Neumann%Neu_Param%neu_index)-1
                index_neu = Tdomain%Neumann%Neu_Param%neu_index(i_neu+1)-1
                write(schar,*) index_neu
                surfacename = "Neumann"//adjustl(schar(1:len_trim(schar)))
                Tdomain%Neumann%NeuSurface(i_neu)%name=Trim(surfacename)
                write (*,1002) surfacename
 
                do i_surf=lbound(Tdomain%sSurfaces,1),ubound(Tdomain%sSurfaces,1)
                   if (trim(Tdomain%sSurfaces(i_surf)%name) == trim(surfacename)) then
 
                     select case (Tdomain%sSurfaces(i_surf)%domain)
                         case (DM_SOLID)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_sl, Tdomain, &
                                            "f", itemp2, itemp2b, itemp2c, itemp)
                         case (DM_FLUID)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fl, Tdomain, &
                                            "f", itemp2, itemp2b, itemp2c, itemp)
                         case (DM_SOLID_PML)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_spml, Tdomain, &
                                            "f", itemp2, itemp2b, itemp2c, itemp)
                         case (DM_FLUID_PML)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fpml, Tdomain, &
                                            "f", itemp2, itemp2b, itemp2c, itemp)
                         case default
                             write (*,*)
                             stop " Unkown material type while reading Neumann surface "
                     end select
 
                     Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_faces = size(itemp2,2)
                     allocate(Tdomain%Neumann%NeuSurface(i_neu)%Neu_face(0:Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_faces-1))
                     do i = 0, Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_faces-1
                         Tdomain%Neumann%NeuSurface(i_neu)%Neu_face(i)%Near_Edges(0:3) = itemp2(:,i+1)
                         Tdomain%Neumann%NeuSurface(i_neu)%Neu_face(i)%Near_Edges_Orient(0:3) = itemp2b(:,i+1)
                     end do
 
                     do i = 0, Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_faces-1
                         Tdomain%Neumann%NeuSurface(i_neu)%Neu_face(i)%Near_Vertices(0:3) = itemp2c(:,i+1)
 
                         Tdomain%Neumann%NeuSurface(i_neu)%Neu_face(i)%Face = itemp(0,i+1)
                     end do
                     deallocate(itemp2, itemp2b, itemp, itemp2c)
                     ! Neumann edges
 
                     select case (Tdomain%sSurfaces(i_surf)%domain)
                          case (DM_SOLID)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_sl, Tdomain, &
                                            "e", itemp, itemp2, itemp2b, itemp2c)
                          case (DM_FLUID)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fl, Tdomain, &
                                            "e", itemp, itemp2, itemp2b, itemp2c)
                          case (DM_SOLID_PML)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_spml, Tdomain, &
                                            "e", itemp, itemp2, itemp2b, itemp2c)
                          case (DM_FLUID_PML)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fpml, Tdomain, &
                                            "e", itemp, itemp2, itemp2b, itemp2c)
                     end select
 
                     Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_edges = size(itemp,2)
                     allocate(Tdomain%Neumann%NeuSurface(i_neu)%Neu_edge(0:Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_edges-1))
                     do i = 0, Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_edges-1
                         Tdomain%Neumann%NeuSurface(i_neu)%Neu_edge(i)%Edge = itemp(0,i+1)
                     end do
                     deallocate(itemp)
                     ! Neumann vertices
 
                     select case (Tdomain%sSurfaces(i_surf)%domain)
                          case (DM_SOLID)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_sl, Tdomain, &
                                            "v", itemp, itemp2, itemp2b, itemp2c)
                          case (DM_FLUID)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fl, Tdomain, &
                                            "v", itemp, itemp2, itemp2b, itemp2c)
                          case (DM_SOLID_PML)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_spml,Tdomain, &
                                            "v", itemp, itemp2, itemp2b, itemp2c)
                          case (DM_FLUID_PML)
                             call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fpml, Tdomain, &
                                            "v", itemp, itemp2, itemp2b, itemp2c)
                     end select
 
                     Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_vertices = size(itemp,2)
                     allocate(Tdomain%Neumann%NeuSurface(i_neu)%Neu_vertex(0:Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_vertices-1))
                     do i = 0, Tdomain%Neumann%NeuSurface(i_neu)%Neu_n_vertices-1
                         Tdomain%Neumann%NeuSurface(i_neu)%Neu_vertex(i)%Vertex = itemp(0,i+1)
                     end do
                     deallocate(itemp)
                  endif
               enddo
             enddo
             
             include 'formats.in'

        end subroutine read_neu_surface
        !!
        !!
        !!
        subroutine read_planeW_surface(Tdomain)
        
        implicit none
        !
        type(domain), intent(inout)          :: Tdomain  
        !
        integer(HID_T)                       :: fid, proc_id, surf_id
        integer                              :: hdferr, ierr, i_planeW, i_surf, i
        integer, allocatable, dimension(:,:) :: itemp2, itemp2b, itemp2c
        integer, allocatable, dimension(:,:) :: itemp
        character(len=12)                    ::  schar
        character(len=100)                   :: surfacename

           do i_planeW=0,0 !Tdomain%n_planeW-1
               write(schar,*) i_planeW
               surfacename = "PlaneWave"//adjustl(schar(1:len_trim(schar)))
               write (*,1002) surfacename

               do i_surf=lbound(Tdomain%sSurfaces,1),ubound(Tdomain%sSurfaces,1)
                  if (trim(Tdomain%sSurfaces(i_surf)%name) == trim(surfacename)) then

                    select case (Tdomain%sSurfaces(i_surf)%domain)
                        case (DM_SOLID)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_sl, Tdomain, &
                                           "f", itemp2, itemp2b, itemp2c, itemp)
                        case (DM_FLUID)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fl, Tdomain, &
                                           "f", itemp2, itemp2b, itemp2c, itemp)
                        case (DM_SOLID_PML)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_spml, Tdomain, &
                                           "f", itemp2, itemp2b, itemp2c, itemp)
                        case (DM_FLUID_PML)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fpml, Tdomain, &
                                           "f", itemp2, itemp2b, itemp2c, itemp)
                        case default
                            write (*,*)
                            stop " Unkown material type while reading plane wave surface "
                    end select

                    Tdomain%sPlaneW%n_faces = size(itemp2,2)
                    allocate(Tdomain%sPlaneW%pface(0:Tdomain%sPlaneW%n_faces-1))
                    do i = 0, Tdomain%sPlaneW%n_faces-1
                        Tdomain%sPlaneW%pface(i)%Near_Edges(0:3) = itemp2(:,i+1)
                        Tdomain%sPlaneW%pface(i)%Orient_Edges(0:3) = itemp2b(:,i+1)
                    end do

                    do i = 0, Tdomain%sPlaneW%n_faces-1
                        Tdomain%sPlaneW%pface(i)%Near_Vertices(0:3) = itemp2c(:,i+1)

                        Tdomain%sPlaneW%pface(i)%Face = itemp(0,i+1)
                    end do
                    deallocate(itemp2, itemp2b, itemp, itemp2c)
                    ! Neumann edges

                    select case (Tdomain%sSurfaces(i_surf)%domain)
                         case (DM_SOLID)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_sl, Tdomain, &
                                           "e", itemp, itemp2, itemp2b, itemp2c)
                         case (DM_FLUID)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fl, Tdomain, &
                                           "e", itemp, itemp2, itemp2b, itemp2c)
                         case (DM_SOLID_PML)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_spml, Tdomain, &
                                           "e", itemp, itemp2, itemp2b, itemp2c)
                         case (DM_FLUID_PML)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fpml, Tdomain, &
                                           "e", itemp, itemp2, itemp2b, itemp2c)
                    end select

                    Tdomain%sPlaneW%n_edges = size(itemp,2)
                    allocate(Tdomain%sPlaneW%pEdge(0:Tdomain%sPlaneW%n_edges-1))
                    do i = 0, Tdomain%sPlaneW%n_edges-1
                        Tdomain%sPlaneW%pEdge(i)%Edge = itemp(0,i+1)
                    end do
                    deallocate(itemp)
                    ! Neumann vertices

                    select case (Tdomain%sSurfaces(i_surf)%domain)
                         case (DM_SOLID)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_sl, Tdomain, &
                                           "v", itemp, itemp2, itemp2b, itemp2c)
                         case (DM_FLUID)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fl, Tdomain, &
                                           "v", itemp, itemp2, itemp2b, itemp2c)
                         case (DM_SOLID_PML)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_spml,Tdomain, &
                                           "v", itemp, itemp2, itemp2b, itemp2c)
                         case (DM_FLUID_PML)
                            call MTExtract(Tdomain%sSurfaces(i_surf)%surf_fpml, Tdomain, &
                                           "v", itemp, itemp2, itemp2b, itemp2c)
                    end select

                    Tdomain%sPlaneW%n_vertices = size(itemp,2)
                    allocate(Tdomain%sPlaneW%pvertex(0:Tdomain%sPlaneW%n_vertices-1))
                    do i = 0, Tdomain%sPlaneW%n_vertices-1
                        Tdomain%sPlaneW%pvertex(i)%Vertex = itemp(0,i+1)
                    end do
                    deallocate(itemp)
                 endif
              enddo
            enddo

            include 'formats.in'

       end subroutine read_planeW_surface
       !!
       !!
       !!
       subroutine MTExtract(Tsurface,Tdomain,type,temp1,temp2,temp3, temp4)
 
         implicit none
         type(domain), intent(in)           :: Tdomain
         type(surf_num), intent(in)         :: Tsurface
         integer, allocatable, intent(inout):: temp1(:,:), temp2(:,:), temp3(:,:), temp4(:,:)
         character, intent(in)              :: type
         integer                            :: i, k, j, n, faces(0:3), edges(0:1), iedg
 
        select case (type)
        case ("f")
           allocate(temp1(0:3,1:Tsurface%n_faces))
           allocate(temp2(0:3,1:Tsurface%n_faces))
           allocate(temp3(0:3,1:Tsurface%n_faces))
           allocate(temp4(0:0,1:Tsurface%n_faces))
           temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0
           do i=1,Tsurface%n_faces
              faces = Tdomain%sFace(Tsurface%if_faces(i-1))%inodes
              iedg = 0; k=0
              block : &
              do
                 edges = Tdomain%sEdge(Tsurface%if_edges(k))%inodes; n = 0
                 do j=0,3
                    if ((faces(j).eq.edges(0)).or.(faces(j).eq.edges(1))) n=n+1
                 enddo
 
                 if (n.eq.2) then
                     !all 4 edges associated to the face
                     temp1(iedg,i)=Tsurface%if_edges(k)
                     !all 4 edges oritation
                     temp2(iedg,i)=Tsurface%if_edge_norm(k)
                     iedg=iedg+1
                  endif
                  k=k+1
                  if ((iedg.ge.4).or.(k.ge.Tsurface%n_edges)) exit block
              enddo block
              !all 4 vertices
              temp3(0:3,i)= Tdomain%sFace(Tsurface%if_faces(i-1))%inodes
              !face id
              temp4(0,i)  = Tsurface%if_faces(i-1)
           enddo
        case ("v")
           allocate(temp1(0:0,1:Tsurface%n_vertices))
           do i=1, Tsurface%n_vertices
              temp1(0,i)=Tsurface%if_vertices(i-1)
           enddo
        case ("e")
            allocate(temp1(0:0,1:Tsurface%n_edges))
            do i=1, Tsurface%n_edges
               temp1(0,i) = Tsurface%if_edges(i-1)
            enddo
     end select
 
     include 'formats.in'

     end subroutine MTExtract 

end module surfSetting
