subroutine shift_to_parameters(rank,n,shift,give_to,take_from,rings)

    implicit none

    integer, intent(in)  :: rank,n,shift
    integer, intent(out) :: give_to,take_from,rings

    give_to = rank+shift
    if(give_to > n-1) give_to = give_to-n
    take_from = rank-shift
    if(take_from < 0) take_from = take_from+n
    if(mod(n,shift) == 0 .and. shift /= 1)then
        rings = shift
    else if(mod(n,n-shift) == 0 .and. shift /= n-1)then
        rings = n-shift
    else if(mod(n,2) == 0 .and. mod(shift,2) == 0)then
        rings = 2
    else
        rings = 1
    endif

end subroutine shift_to_parameters
!------------------------------------------------------------------------
!------------------------------------------------------------------------
subroutine ALGO_COMM(rank,n,rings,ngll_give,ngll_take,give_to,take_from,fact_mult,  &
    TabGiven,TabTaken)
    ! old algo to exchange properties between processors
    use mpi
    implicit none
    integer, intent(in) :: rings,rank,n,ngll_give,ngll_take,give_to,take_from,fact_mult
    real, dimension(0:fact_mult*ngll_give-1), intent(in) :: TabGiven
    real, dimension(0:fact_mult*ngll_take-1), intent(out) :: TabTaken
    integer  :: i,j,code,etiquette,statut

    etiquette = 100

    do i = 0,rings-1
        if(rank == i)then
            if(ngll_give > 0)then
                call MPI_SEND(TabGiven,fact_mult*ngll_give,MPI_DOUBLE_PRECISION,  &
                    give_to,etiquette,MPI_COMM_WORLD,code)
            endif
            if(ngll_take > 0)then
                call MPI_RECV(TabTaken,fact_mult*ngll_take,MPI_DOUBLE_PRECISION,  &
                    take_from,etiquette,MPI_COMM_WORLD,statut,code)
            endif
        else
            do j = 0,n/rings-1
                if(rank == i + j*rings)then
                    if(ngll_take > 0)then
                        call MPI_RECV(TabTaken,fact_mult*ngll_take,MPI_DOUBLE_PRECISION, &
                            take_from,etiquette,MPI_COMM_WORLD,statut,code)
                    endif
                    if(ngll_give > 0)then
                        call MPI_SEND(TabGiven,fact_mult*ngll_give,MPI_DOUBLE_PRECISION, &
                            give_to,etiquette,MPI_COMM_WORLD,code)
                    endif
                endif
            enddo
        endif
    enddo

    return

end subroutine ALGO_COMM
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
