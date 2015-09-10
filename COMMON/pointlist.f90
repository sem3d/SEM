!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file sem_hdf5.f90
!!\brief A data structure to allow unique creation of 3D points
!!\author L. A.
!!\version 1.0
!!\date 2011.11
!!
!! On maintient une hashtable contenant des points 3D. L'ajout d'un point
!! Renvoie un index d'un point existant ou d'un nouveau point
!! Une fois la table construite on peut recalculer un tableau N*(x,y,z)
!!
!<



module point_list

    implicit none
    type :: PointEntry
       real :: x,y,z
       integer :: hash
       integer :: id
       type(PointEntry), pointer :: next
    end type PointEntry

    type :: PointTable
       integer :: n_bins
       real :: scale
       type(PointEntry), pointer, dimension(:) :: entries
       integer :: count
    end type PointTable
contains

    ! Calcul la fonction de hash d'un point 3D.
    !
    ! On fait un simple xor bits a bit des parties entieres des
    ! coordonnees multipliees par scale
    function hash_point(table, x,y,z)
        type(PointTable), intent(in) :: table
        real, intent(in) :: x,y,z
        integer :: hash_point
        !
        integer :: ix, iy, iz, hash
        
        ix = int(x*table%scale)
        iy = int(y*table%scale)
        iz = int(z*table%scale)
        hash = ieor(ix, iy)
        hash = ieor(hash, iz)
        hash_point = hash
    end function hash_point
    
    ! Calcule le numero d'entree dans la table de hash en fonction
    ! du hash et de n_bins. Attention hash peut etre negatif
    function hash_to_index(table, hash)
        type(PointTable), intent(in) :: table
        integer, intent(in) :: hash
        integer :: hash_to_index
        
        hash_to_index = mod(hash, table%n_bins)
        do while (hash_to_index<0)
            hash_to_index = hash_to_index + table%n_bins
        end do
    end function hash_to_index

    subroutine init_point_table(table, nbins)
        type(PointTable), intent(inout) :: table
        integer, intent(in) :: nbins
        !
        integer :: i
        table%n_bins = nbins
        table%count = 0
        table%scale = 123456
        allocate(table%entries(0:nbins-1))
        do i=0,nbins-1
            table%entries(i)%id = -1
            nullify(table%entries(i)%next)
        end do
    end subroutine init_point_table

    subroutine add_point(table, x, y, z, id)
        type(PointTable), intent(inout) :: table
        real, intent(in) :: x,y,z
        integer, intent(out) :: id
        !
        integer :: hash, index
        type(PointEntry), pointer :: en
        !
        hash = hash_point(table, x, y, z)
        index = hash_to_index(table, hash)
        en => table%entries(index)
        do while(en%id /= -1)
            if (en%x==x .and. en%y==y .and. en%z == z) then
                ! Found a matching element
                exit
            end if
            if (.not. associated(en%next)) then
                allocate(en%next)
                en%next%id = -1
                nullify(en%next%next)
            end if
            en => en%next
        end do
        if (en%id==-1) then
            en%id = table%count
            table%count = table%count+1
            en%hash = hash
            en%x = x
            en%y = y
            en%z = z
        end if
        id = en%id
    end subroutine add_point


    subroutine fill_table_array(table, array)
        type(PointTable), intent(in) :: table
        real, dimension(:,:), allocatable, intent(out) :: array
        !
        integer :: i, empty
        type(PointEntry), pointer :: en
        !
        allocate(array(0:2,0:table%count-1))
        !
        empty = 0
        !
        do i=0,table%n_bins-1
            !write(*,*) "Bin:", i
            en => table%entries(i)
            do while(associated(en))
                !write(*,*) "id=", en%id
                if (en%id==-1) then
                    empty = empty + 1
                    exit
                end if
                array(0,en%id) = en%x
                array(1,en%id) = en%y
                array(2,en%id) = en%z
                en => en%next
            end do
        end do
        !DBG
        write(*,*) "Fill rate:", table%n_bins-empty, "/", table%n_bins
        write(*,*) "Avg non-empty bin size:", real(table%count)/real(table%n_bins-empty)
    end subroutine fill_table_array
end module point_list

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
