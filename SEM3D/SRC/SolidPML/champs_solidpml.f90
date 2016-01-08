!>
!!\file champs_solidpml.f90
!!\brief Contient la définition du type champs pour un domaine solide PML
!!
!<

module champs_solidpml

    implicit none

    type :: champssolidpml

        !! Solide PML
        real, dimension(:,:,:), allocatable :: ForcesPML
        real, dimension(:,:,:), allocatable :: VelocPML
        real, dimension(:,:,:), allocatable :: DumpV

    end type champssolidpml

    contains

end module champs_solidpml

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
