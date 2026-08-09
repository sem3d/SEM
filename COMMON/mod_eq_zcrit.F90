!! Módulo compartilhado para o Newmark Modificado
!! Tabela de z_crit pré-calculada para m <= 20 e bissecção para m > 20
module m_mod_eq_zcrit
    use constants
    implicit none

contains

    !> Fator S(z, m) da Equação Modificada Taylor
    function ModifiedEquationS(z, m) result(s)
        implicit none
        real(fpp), intent(in) :: z
        integer, intent(in) :: m
        real(fpp) :: s, term
        integer :: k
        s = -2._fpp + z
        term = z
        do k = 2, m
            term = -term * z / real((2*k)*(2*k - 1), fpp)
            s = s + 2._fpp * term
        end do
    end function ModifiedEquationS

    !> Função de violação do fator de amplificação: max(S, -4 - S)
    function mod_eq_viol(z, m) result(v)
        implicit none
        real(fpp), intent(in) :: z
        integer, intent(in) :: m
        real(fpp) :: v, s
        s = ModifiedEquationS(z, m)
        v = max(s, -4._fpp - s)
    end function mod_eq_viol

    !> Cálculo dinâmico do z_crit via bissecção para ordem m
    function ModifiedEquationZCrit(m) result(zc)
        implicit none
        integer, intent(in) :: m
        real(fpp) :: zc, zlo, zhi, zmid, dz
        integer :: it
        ! Faixa de busca baseada em z_crit = (m+1)^2 aproximadamente
        zhi = (real(m + 2, fpp))**2
        zlo = 0._fpp
        dz = 0.1_fpp
        do while (mod_eq_viol(zhi, m) <= 0._fpp .and. zhi < 10000._fpp)
            zhi = zhi + 10._fpp
        end do
        do it = 1, 100
            zmid = 0.5_fpp * (zlo + zhi)
            if (mod_eq_viol(zmid, m) <= 0._fpp) then
                zlo = zmid
            else
                zhi = zmid
            end if
        end do
        zc = 0.5_fpp * (zlo + zhi)
    end function ModifiedEquationZCrit

    !> Retorna z_crit para a ordem m (tabela direta para m <= 20)
    function get_zcrit(m) result(zc)
        implicit none
        integer, intent(in) :: m
        real(fpp) :: zc

        select case (m)
        case (1)
            zc = 4.00000000000000_fpp
        case (2)
            zc = 12.0000000000000_fpp
        case (4)
            zc = 21.4820847986060_fpp
        case (6)
            zc = 31.0631248041539_fpp
        case (8)
            zc = 40.7061730030588_fpp
        case (10)
            zc = 50.3789421868351_fpp
        case (12)
            zc = 60.0682522774843_fpp
        case (14)
            zc = 69.7674291880921_fpp
        case (16)
            zc = 79.4731818296317_fpp
        case (18)
            zc = 89.1837562810842_fpp
        case (20)
            zc = 98.8979310897711_fpp
        case default
            zc = ModifiedEquationZCrit(m)
        end select
    end function get_zcrit

end module m_mod_eq_zcrit
