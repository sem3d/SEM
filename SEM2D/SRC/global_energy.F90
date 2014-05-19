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
    implicit none
contains

    subroutine global_energy (Tdomain)

        implicit none
        type (domain), intent (INOUT) :: Tdomain

        ! local variables
        integer :: n, mat
        real    :: E_tot, E_k, E_el, Dt

        E_tot = 0.
        E_k  = 0.
        E_el = 0.

        do n=0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%PML) then
                mat = Tdomain%specel(n)%mat_index
                Dt = Tdomain%sSubDomain(mat)%dt
                call compute_Kinetic_Energy (Tdomain%specel(n), Dt, E_k)
                E_tot = E_tot + E_k
                ! Preparation of element displ for computing Elastic Energy
                call Prediction_Elem_Veloc (Tdomain%specel(n))
            endif
        enddo

        ! Computation of the Kinetic Energy for the faces
        do n=0,Tdomain%n_face-1
            if (.not. Tdomain%sFace(n)%PML) then
                mat = Tdomain%sFace(n)%mat_index
                Dt = Tdomain%sSubDomain(mat)%dt
                call compute_Kinetic_Energy_F (Tdomain%sFace(n), Dt, E_k)
                E_tot = E_tot + E_k
                ! Preparation of face displ for computing Elastic Energy
                call Prediction_Face_Veloc (Tdomain%sFace(n))
            endif
        enddo

        ! Computation of the Kinetic Energy for the Vertices
        do n=0,Tdomain%n_vertex-1
            if (.not. Tdomain%sVertex(n)%PML) then
                mat = Tdomain%sVertex(n)%mat_index
                Dt = Tdomain%sSubDomain(mat)%dt
                call compute_Kinetic_Energy_V (Tdomain%sVertex(n), Dt, E_k)
                E_tot = E_tot + E_k
                ! Preparation of vertex displ for computing Elastic Energy
                call Prediction_Vertex_Veloc (Tdomain%sVertex(n))
            endif
        enddo

        do n=0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            ! Computation of the Elastic Energy of the element
            if (.not. Tdomain%specel(n)%PML) then
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

