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

    subroutine global_energy (Tdomain,timelocal)

        implicit none
        type (domain), intent (IN) :: Tdomain
        real, intent (IN)          :: timelocal

        ! local variables
        integer :: n, i, j, mat
        real    :: E_tot, E_k, E_el 

        E_tot = 0.

        do n=0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            ! Computation of the Elastic Energy of the element
            call compute_Elastic_Energy(Tdomain%specel(n), Tdomain%sSubDomain(mat)%hTprime, &
                                        Tdomain%sSubDomain(mat)%hprimez, E_el)
            ! Computation of the Kinetic Energy for the inner nodes of the element
            call compute_Kinetic_Energy (Tdomain%specel(n), E_k)
            E_tot = E_tot + E_el + E_k
        enddo

        ! Computation of the Kinetic Energy for the faces
        do n=0,Tdomain%n_face-1
            E_k = 0.5* sum(1./Tdomain%sFace(n)%MassMat(:)*(Tdomain%sFace(n)%Veloc(:,0)**2 &
                                                          +Tdomain%sFace(n)%Veloc(:,1)**2))
            E_tot = E_tot + E_k
        enddo

        ! Computation of the Kinetic Energy for the Vertices
        do n=0,Tdomain%n_vertex-1
            E_k = 0.5* 1./Tdomain%sVertex(n)%MassMat*( Tdomain%sVertex(n)%Veloc(0)**2 &
                                                      +Tdomain%sVertex(n)%Veloc(1)**2)
            E_tot = E_tot + E_k
        enddo
        
        return
    end subroutine global_energy

end module sglobal_energy
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

