!>
!!\file .F90
!!\brief Algorithme de global_energy
!!\version 1.0
!!\date 15/02/2014
!! La routine global_energy assure le calcul de l'energie globale du systeme :
!! energie_totale = energie_cinetique + energie_deformation_elastique
!<

module sglobal_energy
    use sdomain
    use mpi
    use constants
    implicit none
contains

    subroutine global_energy_generalized (Tdomain)

        implicit none
        type (domain), intent (INOUT) :: Tdomain

        ! local variables
        integer :: n, i, j, mat, type_DG
        real    :: E_tot, E_k, E_el

        E_tot = 0.
        E_k  = 0.
        E_el = 0.
        do n = 0, Tdomain%n_elem-1
            type_DG = Tdomain%specel(n)%Type_DG
            mat = Tdomain%specel(n)%mat_index
            if (type_DG .EQ. GALERKIN_CONT) then
                call get_Displ_fv2el (Tdomain,n)
                call compute_Elastic_Energy (Tdomain%specel(n), Tdomain%sSubDomain(mat)%hTprimex, &
                                        Tdomain%sSubDomain(mat)%hprimez, E_el)
                call compute_Kinetic_Energy (Tdomain%specel(n), E_k)
            else ! Discontinuous Galerkin cases
                call compute_Elastic_Energy_DG (Tdomain%specel(n), E_el)
                call compute_Kinetic_Energy_DG (Tdomain%specel(n), E_k)
            endif
            E_tot = E_tot + E_k + E_el
        enddo

        if (Tdomain%TimeD%rtime == 0) then
            open (51,file = "Total_Energy",status="replace",form="formatted")
        endif
        write(51,*) Tdomain%TimeD%rtime, E_tot

        return
    end subroutine global_energy_generalized



    subroutine global_energy (Tdomain)

        implicit none
        type (domain), intent (INOUT) :: Tdomain

        ! local variables
        integer :: n, i, j, mat
        real    :: E_tot, E_k, E_el

        E_tot = 0.
        E_k  = 0.
        E_el = 0.

        do n=0,Tdomain%n_elem-1
            ! Computation of the Kinetic Energy for the inner nodes of the element
            call compute_Kinetic_Energy (Tdomain%specel(n), E_k)
            E_tot = E_tot + E_k
            ! Preparation of element displ for computing Elastic Energy
            if (.not. Tdomain%specel(n)%PML) call Prediction_Elem_Veloc (Tdomain%specel(n))
        enddo

        ! Computation of the Kinetic Energy for the faces
        do n=0,Tdomain%n_face-1
            E_k = 0.5* sum(1./Tdomain%sFace(n)%MassMat(:)*(Tdomain%sFace(n)%Veloc(:,0)**2 &
                                                          +Tdomain%sFace(n)%Veloc(:,1)**2))
            E_tot = E_tot + E_k
            ! Preparation of face displ for computing Elastic Energy
            if (.not. Tdomain%sFace(n)%PML)  call Prediction_Face_Veloc (Tdomain%sFace(n))
        enddo

        ! Computation of the Kinetic Energy for the Vertices
        do n=0,Tdomain%n_vertex-1
            E_k = 0.5* 1./Tdomain%sVertex(n)%MassMat*( Tdomain%sVertex(n)%Veloc(0)**2 &
                                                      +Tdomain%sVertex(n)%Veloc(1)**2)
            E_tot = E_tot + E_k
            ! Preparation of vertex displ for computing Elastic Energy
            if (.not. Tdomain%sVertex(n)%PML)  call Prediction_Vertex_Veloc (Tdomain%sVertex(n))
        enddo

        do n=0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            ! Computation of the Elastic Energy of the element
            if (.not. Tdomain%specel(n)%PML ) then
                call get_Displ_fv2el (Tdomain,n)
                call compute_Elastic_Energy(Tdomain%specel(n), Tdomain%sSubDomain(mat)%hTprimex, &
                                        Tdomain%sSubDomain(mat)%hprimez, E_el)
                E_tot = E_tot + E_el
            endif
        enddo

        if (Tdomain%TimeD%rtime == 0) then
            open (51,file = "Total_Energy",status="replace",form="formatted")
        endif
        write(51,*) Tdomain%TimeD%rtime, E_tot

        return
    end subroutine global_energy

end module sglobal_energy
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
