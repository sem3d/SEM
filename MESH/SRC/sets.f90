module sets

    implicit none
    type  :: setlist
       integer   :: nr
       integer, dimension(:), pointer  :: elem
    end type setlist

    type  :: near_entity
       integer  :: elem
       type(near_entity), pointer  :: pt => NULL()
    end type near_entity

    type :: near_node
       type(near_entity), pointer  :: ptr => NULL()
    end type near_node

    type :: near_elem
       type(near_entity), pointer  :: ptr => NULL()
    end type near_elem

    type :: near_proc
       integer, pointer, dimension(:)  :: list
       integer  :: nb
    end type near_proc
contains

    subroutine list_union(L1,L2,L)
        ! L1 and L2 must be sorted
        type(setlist), intent(in)   :: L1, L2
        type(setlist), intent(out)  :: L
        integer  :: i, i1,i2,p


        i = 0 ; i1 = 0; i2 = 0

        do
            if(i1 >= L1%nr .or. i2 >= L2%nr) exit
            p = L1%elem(i1)-L2%elem(i2)
            if(p < 0)then
                L%elem(i) = L1%elem(i1)
                i = i+1 ; i1 = i1+1
            else if(p > 0)then
                L%elem(i) = L2%elem(i2)
                i = i+1 ; i2 = i2+1
            else
                L%elem(i) = L1%elem(i1)
                i = i+1 ; i1 = i1+1 ; i2 = i2+1
            end if
        end do
        if(i1 >= L1%nr)then
            do
                if(i2 >= L2%nr) exit
                L%elem(i) = L2%elem(i2)
                i = i+1 ; i2 = i2+1
            end do

        end if

        if(i2 >= L2%nr)then
            do
                if(i1 >= L1%nr) exit
                L%elem(i) = L1%elem(i1)
                i = i+1 ; i1 = i1+1
            end do
        end if

        L%nr = i

    end subroutine list_union
    !---------------------------------------------
    !------------------------------------------------
    recursive subroutine list_entity(ptrent)
        type(near_entity), pointer  :: ptrent

        if(associated(ptrent%pt)) call list_entity(ptrent%pt)
        print*,"HE",ptrent%elem

    end subroutine list_entity
    !-----------------------------------------------
    logical function is_in_the_list(ptr,ndat)
        implicit none
        type(near_entity), pointer   :: ptr
        integer, intent(in)          :: ndat
        type(near_entity), pointer   :: ptloc

        ptloc => ptr

        do while(associated(ptloc))
            if(ptloc%elem == ndat)then
                is_in_the_list = .true.
                return
            else
                ptloc => ptloc%pt
            end if
        end do

        is_in_the_list = .false.
        return

    end function is_in_the_list
    !-------------------------------------------------
    subroutine entity_sort(ptr)
        implicit none
        type(near_entity), pointer    :: ptr
        type(near_entity), pointer    :: ptloc1, ptloc2, ptrinit
        integer   :: tmp

        ptloc1 => ptr
        ptrinit => ptloc1   ! 1er de la liste

        do while(associated(ptloc1%pt))
            ptloc2 => ptloc1%pt
            do while(associated(ptloc2))
                if(ptloc2%elem < ptloc1%elem)then
                    tmp = ptloc2%elem
                    ptloc2%elem = ptloc1%elem
                    ptloc1%elem = tmp
                end if
                ptloc2 => ptloc2%pt
            end do
            ptloc1 => ptloc1%pt
        end do

        ptr => ptrinit

    end subroutine entity_sort
    !-------------------------------------------------
    subroutine push_entity(ptrent,ndat)
        implicit none
        type(near_entity), pointer   :: ptrent
        integer, intent(in)          :: ndat
        type(near_entity), pointer   :: ptloc

        allocate(ptloc)
        ptloc%elem = ndat
        ptloc%pt => ptrent
        ptrent => ptloc

    end subroutine push_entity
    !------------------------------------------------
    subroutine entity_union(list1,list2)
        ! adds list2 to the first one
        implicit none
        type(near_entity), pointer  :: list1, list2
        type(near_entity), pointer  :: listloc

        if(.not. associated(list1) .and. .not. associated(list2)) return
        listloc => list2

        do while(associated(listloc))
            if(.not. is_in_the_list(list1,listloc%elem))then
                call push_entity(list1,listloc%elem)
            end if
            listloc => listloc%pt
        end do

    end subroutine entity_union
    !------------------------------------------------
    subroutine entity_intersect(list1,list2,list)
        ! intersection of 2 lists that must be sorted
        implicit none
        type(near_entity), pointer   :: list1, list2, list
        type(near_entity), pointer   :: listloc1, listloc2

        list => NULL()
        if(.not. associated(list1) .or. .not. associated(list2)) return
        listloc1 => list1 ; listloc2 => list2

        inter_loop: do while(associated(listloc1))
            if(listloc1%elem > listloc2%elem)then
                listloc2 => listloc2%pt
            else if(listloc1%elem < listloc2%elem)then
                listloc1 => listloc1%pt
            else  !! equality
                call push_entity(list,listloc1%elem)
                listloc1 => listloc1%pt ; listloc2 => listloc2%pt
            end if
            if(.not. associated(listloc2)) exit inter_loop
        end do inter_loop

    end subroutine entity_intersect
    !------------------------------------------------
    subroutine find_near_elem(Ipoint,tab_ptr_node,tab_ptr_elem)
        implicit none
        integer, dimension(0:,0:), intent(in)         :: Ipoint
        type(near_node), dimension(0:), intent(in)    :: tab_ptr_node
        type(near_elem), dimension(0:), intent(inout) :: tab_ptr_elem
        integer  :: n, i, j

        do n = 0, size(Ipoint,2)-1
            do i = 0, size(Ipoint,1)-1
                j = Ipoint(i,n)
                call entity_union(tab_ptr_elem(n)%ptr,tab_ptr_node(j)%ptr)
            end do
        end do

    end subroutine find_near_elem
    !------------------------------------------------
    integer function already_in(ii,tab,nn)
        integer :: ii,nn
        integer, dimension(0:), intent(in) :: tab
        integer :: i
        already_in = -1
        do i = 0,nn-1
            if (ii==tab(i)) already_in = i
        enddo
    end function already_in
    !------------------------------------------------
    !--------------------------------------
    subroutine sort(vector, length)

        implicit none

        integer, intent(IN) :: length
        integer, dimension(0:length-1), intent(INOUT) :: vector
        integer :: n, ok, i, tmp

        n = length-1
        ok = 1
        do while (ok==1)
            ok = 0
            do i = 0,n-1
                if (vector(i)>vector(i+1)) then
                    tmp = vector(i+1)
                    vector(i+1) = vector(i)
                    vector(i) = tmp
                    ok = 1
                endif
            enddo
            n = n-1
        enddo

        return
    end subroutine sort
    !---------------------------------------
    !----------------------------
    subroutine union_sets_elem(setnode,setelem,Ipoint,near_size)


        type(setlist), intent(in)  :: setnode(0:)
        type(setlist), intent(inout) :: setelem(0:)
        integer, intent(in)        :: Ipoint(0:,0:)
        type(setlist)              :: set_inter, set_inter_out
        integer                    :: n,i,nnod, near_size

        allocate(set_inter%elem(0:near_size-1))
        allocate(set_inter_out%elem(0:near_size-1))

        do n = 0, size(setelem)-1
            set_inter_out%elem(0:) = 0 ; set_inter_out%nr = 0
            do i = 0,7
                nnod = Ipoint(i,n)
                set_inter%nr = set_inter_out%nr ; set_inter%elem(0:) = set_inter_out%elem(0:)
                call list_union(setnode(nnod),set_inter,set_inter_out)
            end do
            setelem(n)%nr = set_inter_out%nr ; setelem(n)%elem(0:) = set_inter_out%elem(0:)
        end do

    end subroutine union_sets_elem
    !---------------------------------------
    !----------------------------
    subroutine find_near_proc(n_elem,near_elem_set,part,elem_near_proc)
        implicit none
        integer, intent(in)   ::  n_elem
        integer, intent(in), dimension(0:n_elem-1)  :: part
        type(near_elem), dimension(0:n_elem-1), intent(in) :: near_elem_set
        type(near_proc), dimension(0:n_elem-1), intent(out) :: elem_near_proc
        integer   :: nel,n,proc,i,j,k,this_proc
        integer, dimension(0:49)  :: tr_elem
        type(near_entity), pointer  :: near_neighb => NULL()

        do n = 0,n_elem-1
            i = 0 ; tr_elem = 0 ; this_proc = part(n)
            near_neighb => near_elem_set(n)%ptr  ! pointer on the linked list of nearest elements
            do while(associated(near_neighb))
                proc = part(near_neighb%elem)
                if(proc <= this_proc)then
                    near_neighb => near_neighb%pt
                    cycle
                end if
                k = already_in(proc,tr_elem,i)
                if(k < 0)then
                    tr_elem(i) = proc
                    i = i+1
                end if
                near_neighb => near_neighb%pt
            end do
            elem_near_proc(n)%nb = i
            allocate(elem_near_proc(n)%list(0:i-1))
            do j = 0,i-1
                elem_near_proc(n)%list(j) = tr_elem(j)
            end do
        end do

    end subroutine find_near_proc
    !----------------------------
                                                                       

end module sets
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
