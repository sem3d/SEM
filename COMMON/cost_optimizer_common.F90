!! Módulo compartilhado para ordenação e avaliação de custo espectral
module m_cost_optimizer_common
    use constants
    use m_mod_eq_zcrit
    implicit none

contains

    !> Algoritmo In-Place Quicksort para arrays de números reais real(fpp)
    recursive subroutine quicksort_real(a, first, last)
        implicit none
        real(fpp), dimension(:), intent(inout) :: a
        integer, intent(in) :: first, last
        integer :: i, j
        real(fpp) :: x, temp

        x = a((first + last) / 2)
        i = first
        j = last
        do
            do while (a(i) < x)
                i = i + 1
            end do
            do while (x < a(j))
                j = j - 1
            end do
            if (i >= j) exit
            temp = a(i); a(i) = a(j); a(j) = temp
            i = i + 1
            j = j - 1
        end do
        if (first < i - 1) call quicksort_real(a, first, i - 1)
        if (j + 1 < last)  call quicksort_real(a, j + 1, last)
    end subroutine quicksort_real

    !> Determina as ordens m_base para um candidato dt_target em elementos locais
    subroutine compute_element_base_orders(dt_elem_loc, dt_target_cand, m_base_local, cand_feasible)
        implicit none
        real(fpp), dimension(:), intent(in) :: dt_elem_loc
        real(fpp), intent(in) :: dt_target_cand
        integer, dimension(:), intent(out) :: m_base_local
        logical, intent(out) :: cand_feasible
        !
        integer :: n, n_local, m_candidate
        real(fpp) :: r_req

        n_local = size(dt_elem_loc)
        cand_feasible = .true.

        do n = 1, n_local
            if (dt_elem_loc(n) <= 0._fpp .or. dt_target_cand <= dt_elem_loc(n)) then
                m_base_local(n) = 1
            else
                r_req = 4._fpp * (dt_target_cand / dt_elem_loc(n))**2
                m_candidate = 1
                do while (m_candidate <= 20 .and. get_zcrit(m_candidate) < r_req)
                    if (m_candidate == 1) then
                        m_candidate = 2
                    else
                        m_candidate = m_candidate + 2
                    end if
                end do
                if (m_candidate > 20) then
                    cand_feasible = .false.
                    exit
                end if
                m_base_local(n) = m_candidate
            end if
        end do
    end subroutine compute_element_base_orders

end module m_cost_optimizer_common
