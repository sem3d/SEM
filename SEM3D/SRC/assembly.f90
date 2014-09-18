module assembly
    use orientation

contains

    ! Forces assembly

    subroutine get_Forces_Elem2Edge(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ngllx,nglly,ngllz,ngll,ne,nne,orient_e

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do ne = 0,11
            nne = Tdomain%specel(n)%Near_Edges(ne)
            orient_e = Tdomain%specel(n)%Orient_Edges(ne)
            ngll = Tdomain%sEdge(nne)%ngll
            ! now we call the general assemblage routine
            call get_VectProperty_Elem2Edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%Forces(:,0:2),Tdomain%specel(n)%sl%Forces(:,:,:,0:2))
            if(Tdomain%sEdge(nne)%PML)then
                call get_VectProperty_Elem2Edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                    Tdomain%sEdge(nne)%spml%Forces1(:,0:2),Tdomain%specel(n)%slpml%Forces1(:,:,:,0:2))
                call get_VectProperty_Elem2Edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                    Tdomain%sEdge(nne)%spml%Forces2(:,0:2),Tdomain%specel(n)%slpml%Forces2(:,:,:,0:2))
                call get_VectProperty_Elem2Edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                    Tdomain%sEdge(nne)%spml%Forces3(:,0:2),Tdomain%specel(n)%slpml%Forces3(:,:,:,0:2))
            end if
        end do

        return

    end subroutine get_Forces_Elem2Edge
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    subroutine get_ForcesFl_Elem2Edge(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ngllx,nglly,ngllz,ngll,ne,nne,orient_e

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do ne = 0,11
            nne = Tdomain%specel(n)%Near_Edges(ne)
            orient_e = Tdomain%specel(n)%Orient_Edges(ne)
            ngll = Tdomain%sEdge(nne)%ngll
            ! now we call the general assemblage routine
            call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%ForcesFl(:),Tdomain%specel(n)%fl%ForcesFl(:,:,:))
            if(Tdomain%sEdge(nne)%PML)then
                call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                    Tdomain%sEdge(nne)%spml%ForcesFl1(:),Tdomain%specel(n)%flpml%ForcesFl1(:,:,:))
                call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                    Tdomain%sEdge(nne)%spml%ForcesFl2(:),Tdomain%specel(n)%flpml%ForcesFl2(:,:,:))
                call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                    Tdomain%sEdge(nne)%spml%ForcesFl3(:),Tdomain%specel(n)%flpml%ForcesFl3(:,:,:))
            end if
        end do

        return

    end subroutine get_ForcesFl_Elem2Edge

    subroutine get_Forces_Elem2Face(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ngllx,nglly,ngllz,ngll1,ngll2,nf,nnf,orient_f

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2
            ! now we call the general assemblage routine
            call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%Forces(:,:,0:2),Tdomain%specel(n)%sl%Forces(:,:,:,0:2))
            if(Tdomain%sFace(nnf)%PML)then
                call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                    Tdomain%sFace(nnf)%spml%Forces1(:,:,0:2),Tdomain%specel(n)%slpml%Forces1(:,:,:,0:2))
                call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                    Tdomain%sFace(nnf)%spml%Forces2(:,:,0:2),Tdomain%specel(n)%slpml%Forces2(:,:,:,0:2))
                call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                    Tdomain%sFace(nnf)%spml%Forces3(:,:,0:2),Tdomain%specel(n)%slpml%Forces3(:,:,:,0:2))
            end if
        enddo

        return

    end subroutine get_Forces_Elem2Face
    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    subroutine get_ForcesFl_Elem2Face(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ngllx,nglly,ngllz,ngll1,ngll2,nf,nnf,orient_f

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2
            ! now we call the general assemblage routine
            call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%ForcesFl(:,:),Tdomain%specel(n)%fl%ForcesFl(:,:,:))
            if(Tdomain%sFace(nnf)%PML)then
                call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                    Tdomain%sFace(nnf)%spml%ForcesFl1(:,:),Tdomain%specel(n)%flpml%ForcesFl1(:,:,:))
                call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                    Tdomain%sFace(nnf)%spml%ForcesFl2(:,:),Tdomain%specel(n)%flpml%ForcesFl2(:,:,:))
                call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                    Tdomain%sFace(nnf)%spml%ForcesFl3(:,:),Tdomain%specel(n)%flpml%ForcesFl3(:,:,:))
            end if
        enddo

        return

    end subroutine get_ForcesFl_Elem2Face


    subroutine get_Forces_Elem2Vertex(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ngllx,nglly,ngllz,nv,nnv


        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            ! now we call the general assemblage routine
            call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%Forces(0:2),Tdomain%specel(n)%sl%Forces(:,:,:,0:2))
            if(Tdomain%sVertex(nnv)%PML)then
                call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                    Tdomain%sVertex(nnv)%spml%Forces1(0:2),Tdomain%specel(n)%slpml%Forces1(:,:,:,0:2))
                call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                    Tdomain%sVertex(nnv)%spml%Forces2(0:2),Tdomain%specel(n)%slpml%Forces2(:,:,:,0:2))
                call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                    Tdomain%sVertex(nnv)%spml%Forces3(0:2),Tdomain%specel(n)%slpml%Forces3(:,:,:,0:2))
            end if
        enddo

        return

    end subroutine get_Forces_Elem2Vertex
    !---------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------
    subroutine get_ForcesFl_Elem2Vertex(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ngllx,nglly,ngllz,nv,nnv


        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            ! now we call the general assemblage routine
            call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%ForcesFl,Tdomain%specel(n)%fl%ForcesFl(:,:,:))
            if(Tdomain%sVertex(nnv)%PML)then
                call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                    Tdomain%sVertex(nnv)%spml%ForcesFl1,Tdomain%specel(n)%flpml%ForcesFl1(:,:,:))
                call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                    Tdomain%sVertex(nnv)%spml%ForcesFl2,Tdomain%specel(n)%flpml%ForcesFl2(:,:,:))
                call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                    Tdomain%sVertex(nnv)%spml%ForcesFl3,Tdomain%specel(n)%flpml%ForcesFl3(:,:,:))
            end if
        enddo

        return

    end subroutine get_ForcesFl_Elem2Vertex



    ! Disp assembly
    ! =======================================

    subroutine get_Displ_Edge2Elem(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ne,ngllx,nglly,ngllz,ngll,nne,orient_e


        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do ne = 0,11
            nne = Tdomain%specel(n)%Near_Edges(ne)
            orient_e = Tdomain%specel(n)%Orient_Edges(ne)
            ngll = Tdomain%sEdge(nne)%ngll
            ! now we call the general deassemblage routine
            call get_VectProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%Forces(:,0:2),Tdomain%specel(n)%sl%Forces(:,:,:,0:2))
        end do

        return

    end subroutine get_Displ_Edge2Elem
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    subroutine get_Phi_Edge2Elem(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ne,ngllx,nglly,ngllz,ngll,nne,orient_e


        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do ne = 0,11
            nne = Tdomain%specel(n)%Near_Edges(ne)
            orient_e = Tdomain%specel(n)%Orient_Edges(ne)
            ngll = Tdomain%sEdge(nne)%ngll
            ! now we call the general deassemblage routine
            call get_ScalarProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%ForcesFl(:),Tdomain%specel(n)%fl%ForcesFl(:,:,:))
        end do

        return

    end subroutine get_Phi_Edge2Elem

    subroutine get_Displ_Face2Elem(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: nf,ngllx,nglly,ngllz,ngll1,ngll2,nnf,orient_f

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2
            ! now we call the general deassemblage routine
            call get_VectProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%Forces(:,:,0:2),Tdomain%specel(n)%sl%Forces(:,:,:,0:2))
        enddo

        return

    end subroutine get_Displ_Face2Elem
    !---------------------------------------------------------

    !---------------------------------------------------------
    subroutine get_Phi_Face2Elem(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: nf,ngllx,nglly,ngllz,ngll1,ngll2,nnf,orient_f

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2
            ! now we call the general deassemblage routine
            call get_ScalarProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%ForcesFl(:,:),Tdomain%specel(n)%fl%ForcesFl(:,:,:))
        enddo

        return

    end subroutine get_Phi_Face2Elem

    subroutine get_Displ_Vertex2Elem(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: nv, ngllx,nglly,ngllz,nnv

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            ! now we call the general deassemblage routine
            call get_VectProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%Forces(0:2),Tdomain%specel(n)%sl%Forces(:,:,:,0:2))
        enddo

        return
    end subroutine get_Displ_Vertex2Elem

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    subroutine get_Phi_Vertex2Elem(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: nv, ngllx,nglly,ngllz,nnv

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            ! now we call the general deassemblage routine
            call get_ScalarProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%ForcesFl,Tdomain%specel(n)%fl%ForcesFl(:,:,:))
        enddo

        return
    end subroutine get_Phi_Vertex2Elem


    ! Mass assembly
    ! =============================================

    subroutine get_Mass_Elem2Edge(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ngllx,nglly,ngllz,ngll,ne,nne,orient_e

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do ne = 0,11
            nne = Tdomain%specel(n)%Near_Edges(ne)
            orient_e = Tdomain%specel(n)%Orient_Edges(ne)
            ngll = Tdomain%sEdge(nne)%ngll
            ! now we call the general assemblage routine
            call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%MassMat,Tdomain%specel(n)%MassMat)
            if(Tdomain%sEdge(nne)%PML)then
                call get_VectProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                    Tdomain%sEdge(nne)%spml%DumpMass(:,0:2),Tdomain%specel(n)%xpml%DumpMass(:,:,:,0:2))
                if(Tdomain%sEdge(nne)%FPML)then
                    call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                        Tdomain%sEdge(nne)%spml%Ivx(:),Tdomain%specel(n)%slpml%Ivx(:,:,:))
                    call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                        Tdomain%sEdge(nne)%spml%Ivy(:),Tdomain%specel(n)%slpml%Ivy(:,:,:))
                    call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                        Tdomain%sEdge(nne)%spml%Ivz(:),Tdomain%specel(n)%slpml%Ivz(:,:,:))
                end if
            end if

        enddo

        return

    end subroutine get_Mass_Elem2Edge

    subroutine get_Mass_Elem2Face(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ngllx,nglly,ngllz,ngll1,ngll2,nf,nnf,orient_f

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2

            ! now we call the general assemblage routine
            call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%MassMat(:,:),Tdomain%specel(n)%MassMat(:,:,:))
            if(Tdomain%sFace(nnf)%PML)then
                call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                    Tdomain%sFace(nnf)%spml%DumpMass(:,:,0:2),Tdomain%specel(n)%xpml%DumpMass(:,:,:,0:2))
                if(Tdomain%sFace(nnf)%FPML)then
                    call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                        Tdomain%sFace(nnf)%spml%Ivx(:,:),Tdomain%specel(n)%slpml%Ivx(:,:,:))
                    call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                        Tdomain%sFace(nnf)%spml%Ivy(:,:),Tdomain%specel(n)%slpml%Ivy(:,:,:))
                    call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                        Tdomain%sFace(nnf)%spml%Ivz(:,:),Tdomain%specel(n)%slpml%Ivz(:,:,:))
                end if
            end if

        enddo

        return

    end subroutine get_Mass_Elem2Face


    subroutine get_Mass_Elem2Vertex(Tdomain,n)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        integer :: ngllx,nglly,ngllz,nv,nnv


        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            ! now we call the general assemblage routine
            call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%MassMat,Tdomain%specel(n)%MassMat(:,:,:))
            if(Tdomain%sVertex(nnv)%PML)then
                call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                    Tdomain%sVertex(nnv)%spml%DumpMass(0:2),Tdomain%specel(n)%xpml%DumpMass(:,:,:,0:2))
                if(Tdomain%sVertex(nnv)%FPML)then
                    call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                        Tdomain%sVertex(nnv)%spml%Ivx,Tdomain%specel(n)%slpml%Ivx(:,:,:))
                    call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                        Tdomain%sVertex(nnv)%spml%Ivy,Tdomain%specel(n)%slpml%Ivy(:,:,:))
                    call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                        Tdomain%sVertex(nnv)%spml%Ivz,Tdomain%specel(n)%slpml%Ivz(:,:,:))
                end if
            end if

        enddo

        return

    end subroutine get_Mass_Elem2Vertex



    ! PML Assembly routines
    ! ==============================================


    subroutine get_PMLprediction_f2el(Tdomain,n,bega,dt)

        use sdomain
        implicit none

        type(Domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        real, intent(in) :: dt, bega
        integer :: nf,ngllx,nglly,ngllz,ngll1,ngll2,nnf,orient_f

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2
            ! now we call the general deassemblage routine
            call get_VectProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,   &
                ngll2,Tdomain%sFace(nnf)%Veloc(:,:,0:2)+              &
                dt*(0.5-bega)*Tdomain%sFace(nnf)%Accel(:,:,0:2),           &
                Tdomain%specel(n)%sl%Forces(:,:,:,0:2))
        enddo

        return

    end subroutine get_PMLprediction_f2el
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    subroutine get_PMLprediction_f2el_fl(Tdomain,n,bega,dt)

        use sdomain
        implicit none

        type(Domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        real, intent(in) :: dt, bega
        integer :: nf,ngllx,nglly,ngllz,ngll1,ngll2,nnf,orient_f

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2
            ! now we call the general deassemblage routine
            call get_ScalarProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,       &
                Tdomain%sFace(nnf)%VelPhi(:,:)+dt*(0.5-bega)*Tdomain%sFace(nnf)%AccelPhi(:,:), &
                Tdomain%specel(n)%fl%ForcesFl(:,:,:))
        enddo

        return

    end subroutine get_PMLprediction_f2el_fl

    subroutine get_PMLprediction_v2el (Tdomain, n, bega, dt)

        use sdomain

        implicit none

        type (Domain), intent (INOUT) :: Tdomain
        integer, intent (IN) :: n
        real, intent(IN) :: dt, bega

        integer :: nv, ngllx, nglly, ngllz, nnv

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            ! now we call the general deassemblage routine
            call get_VectProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,                             &
                Tdomain%sVertex(nnv)%Veloc(0:2)+dt*(0.5-bega)*Tdomain%sVertex(nnv)%Accel(0:2),  &
                Tdomain%specel(n)%sl%Forces(:,:,:,0:2))
        enddo

        return
    end subroutine get_PMLprediction_v2el
    !----------------------------------------------------------

    !----------------------------------------------------------
    subroutine get_PMLprediction_v2el_fl(Tdomain,n,bega,dt)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        real, intent(in)   :: bega,dt
        integer :: nv, ngllx,nglly,ngllz,nnv

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            ! now we call the general deassemblage routine
            call get_ScalarProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,                     &
                Tdomain%sVertex(nnv)%VelPhi+dt*(0.5-bega)*Tdomain%sVertex(nnv)%AccelPhi,  &
                Tdomain%specel(n)%fl%ForcesFl(:,:,:))
        enddo

        return
    end subroutine get_PMLprediction_v2el_fl


    subroutine get_PMLprediction_e2el(Tdomain,n,bega,dt)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        real, intent(in)  :: bega,dt
        integer :: ne,ngllx,nglly,ngllz,ngll,nne,orient_e

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do ne = 0,11
            nne = Tdomain%specel(n)%Near_Edges(ne)
            orient_e = Tdomain%specel(n)%Orient_Edges(ne)
            ngll = Tdomain%sEdge(nne)%ngll
            ! now we call the general deassemblage routine
            call get_VectProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,            &
                Tdomain%sEdge(nne)%Veloc(:,0:2)+dt*(0.5-bega)*Tdomain%sEdge(nne)%Accel(:,0:2), &
                Tdomain%specel(n)%sl%Forces(:,:,:,0:2))
        end do

        return

    end subroutine get_PMLprediction_e2el
    !-----------------------------------------------------------------
    !-----------------------------------------------------------------
    subroutine get_PMLprediction_e2el_fl(Tdomain,n,bega,dt)

        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n
        real, intent(in)  :: bega,dt
        integer :: ne,ngllx,nglly,ngllz,ngll,nne,orient_e

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do ne = 0,11
            nne = Tdomain%specel(n)%Near_Edges(ne)
            orient_e = Tdomain%specel(n)%Orient_Edges(ne)
            ngll = Tdomain%sEdge(nne)%ngll
            ! now we call the general deassemblage routine
            call get_ScalarProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,          &
                Tdomain%sEdge(nne)%VelPhi(:)+dt*(0.5-bega)*Tdomain%sEdge(nne)%AccelPhi(:), &
                Tdomain%specel(n)%fl%ForcesFl(:,:,:))
        end do

        return

    end subroutine get_PMLprediction_e2el_fl

end module assembly
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
