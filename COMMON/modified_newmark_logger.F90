!! Módulo compartilhado para saída de logs e estatísticas do Newmark Modificado
module m_modified_newmark_logger
    use constants
    implicit none

contains

    !> Exibe cabeçalho estilizado da otimização de custo
    subroutine log_header_modified_newmark(dim_str)
        implicit none
        character(len=*), intent(in) :: dim_str
        write(*, fmt='(A)') ""
        write(*, fmt='(A)') "================================================================================"
        write(*, fmt='(A)') "      NEWMARK MODIFICADO REGIONAL (" // trim(dim_str) // ") - OTIMIZAÇÃO DE CUSTO & BFS"
        write(*, fmt='(A)') "================================================================================"
    end subroutine log_header_modified_newmark

    !> Exibe linha da avaliação de candidato a dt_target
    subroutine log_candidate_eval(c, dt_target, weight, cost, is_optimal, is_feasible)
        implicit none
        integer, intent(in) :: c
        real(fpp), intent(in) :: dt_target, weight, cost
        logical, intent(in) :: is_optimal, is_feasible

        if (.not. is_feasible) then
            write(*, fmt='(A,I2,A,1PE12.5,A)') &
                " Candidate #", c, " | dt_target = ", dt_target, &
                " | INVIÁVEL (exige ordem m > 20)"
        else if (is_optimal) then
            write(*, fmt='(A,I2,A,1PE12.5,A,F10.1,A,1PE12.5,A)') &
                " Candidate #", c, " | dt_target = ", dt_target, &
                " | Peso Malha = ", weight, " | Custo = ", cost, " [MELHOR]"
        else
            write(*, fmt='(A,I2,A,1PE12.5,A,F10.1,A,1PE12.5)') &
                " Candidate #", c, " | dt_target = ", dt_target, &
                " | Peso Malha = ", weight, " | Custo = ", cost
        end if
    end subroutine log_candidate_eval

    !> Exibe relatório da configuração ótima selecionada
    subroutine log_optimal_selection(dt_target_opt, min_cost, cost_base, speedup_est)
        implicit none
        real(fpp), intent(in) :: dt_target_opt, min_cost, cost_base, speedup_est

        write(*, fmt='(A)') ""
        write(*, fmt='(A)') " ================================================================================"
        write(*, fmt='(A)') " CONFIGURAÇÃO ÓTIMA SELECIONADA:"
        write(*, fmt='(A,1PE12.5,A)') "  --> Passo de Tempo Global Alvo (dt_alvo) = ", dt_target_opt, " s"
        write(*, fmt='(A,1PE12.5,A,1PE12.5,A)') "  --> Custo Computacional Estimado        = ", min_cost, " (vs Base: ", cost_base, ")"
        write(*, fmt='(A,F6.2,A)') "  --> Speedup Computacional Estimado      = ", speedup_est, "x"
        write(*, fmt='(A)') " ================================================================================"
        write(*, fmt='(A)') ""
    end subroutine log_optimal_selection

    !> Exibe o histograma da distribuição de ordens espectrais m
    subroutine log_order_histogram(hist_g, n_total_g, max_order)
        implicit none
        integer, dimension(0:), intent(in) :: hist_g
        integer, intent(in) :: n_total_g, max_order
        integer :: m
        real(fpp) :: pct

        write(*, fmt='(A)') " --- Distribuição Final de Ordens Espectrais m na Malha ---"
        do m = 1, min(max_order, 20)
            if (m == 1 .or. hist_g(m) > 0) then
                pct = 100._fpp * (real(hist_g(m), fpp) / max(real(n_total_g, fpp), 1._fpp))
                write(*, fmt='(A,I2,A,I8,A,F5.1,A)') &
                    "  [Ordem m = ", m, "] : ", hist_g(m), " elementos (", pct, " %)"
            end if
        end do
        write(*, fmt='(A)') ""
    end subroutine log_order_histogram

end module m_modified_newmark_logger
