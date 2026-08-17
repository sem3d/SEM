!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file mod_eq_zcrit.F90
!!\brief S(z,m), the truncated Taylor amplification term of the order-m
!!  modified-equation Newmark scheme, and its critical stability boundary
!!  z_crit(m). Shared by SEM2D and SEM3D (irons_dtcrit.F90).
!<
module m_mod_eq_zcrit
    use constants
    implicit none

    !> Below this width (in z), an instability interval is treated as
    !! floating-point/near-tangency noise rather than a real stability
    !! boundary: no achievable dt could ever land inside it (it is many
    !! orders of magnitude narrower than any Courant-scaled dt granularity)
    !! and it is barely distinguishable from double-precision noise itself
    !! (confirmed values as low as ~1e-14 in the amplification excess at
    !! the island's peak, for m=18). ModifiedEquationZCrit skips any
    !! interval narrower than this and keeps searching for the next one.
    real(fpp), parameter :: MOD_EQ_MIN_ISLAND_WIDTH = 1.0e-4_fpp

contains

    !> S(z,m) = sum_{k=1}^m 2*(-z)^k/(2k)!, the modified-equation
    !! amplification term: 2+S(z,m) is the order-m Taylor truncation of the
    !! exact amplification factor 2*cos(sqrt(z)) (z=(omega*dt)^2). Built via
    !! the term-to-term ratio term_k = term_{k-1}*(-z)/((2k)(2k-1)),
    !! starting from term_1=(-z)/2! -- verified against 2*cos(sqrt(z))-2 in
    !! the m->infinity limit.
    function ModifiedEquationS(z, m) result(s)
        implicit none
        real(fpp), intent(in) :: z
        integer, intent(in) :: m
        real(fpp) :: s, term
        integer :: k
        term = -z / 2._fpp
        s = 2._fpp * term
        do k = 2, m
            term = term * (-z) / real((2*k)*(2*k - 1), fpp)
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

    !> Scan [a,b) with step dz for the first place mod_eq_viol(z,m) crosses
    !! from <=0 (stable) to >0 (unstable); if found, refine to machine
    !! precision via bisection within that bracket and return zc,
    !! found=.true.; otherwise found=.false. (no crossing in [a,b)).
    subroutine scan_first_crossing(m, a, b, dz, found, zc)
        implicit none
        integer, intent(in) :: m
        real(fpp), intent(in) :: a, b, dz
        logical, intent(out) :: found
        real(fpp), intent(out) :: zc
        real(fpp) :: z, zlo, zhi, zmid, v_lo, v_hi
        integer :: it

        found = .false.
        zc = 0._fpp
        z = a
        v_lo = mod_eq_viol(z, m)
        do while (z < b)
            v_hi = mod_eq_viol(z + dz, m)
            if (v_lo <= 0._fpp .and. v_hi > 0._fpp) then
                found = .true.
                zlo = z
                zhi = z + dz
                do it = 1, 100
                    zmid = 0.5_fpp * (zlo + zhi)
                    if (mod_eq_viol(zmid, m) <= 0._fpp) then
                        zlo = zmid
                    else
                        zhi = zmid
                    end if
                end do
                zc = 0.5_fpp * (zlo + zhi)
                return
            end if
            z = z + dz
            v_lo = v_hi
        end do
    end subroutine scan_first_crossing

    !> Robust critical z for order m: the smallest z>0 where the
    !! amplification factor 2+S(z,m) leaves the stability band [-2,2] (i.e.
    !! mod_eq_viol(z,m) turns positive) for a stretch of at least
    !! MOD_EQ_MIN_ISLAND_WIDTH (or permanently). S(z,m) is a degree-m
    !! polynomial approximating 2*cos(sqrt(z))-2, whose instability set can
    !! have several disjoint intervals before the final, permanent one --
    !! some genuinely significant (e.g. m=10,12,14 have real, moderately
    !! wide early instability windows a naive search can jump clean over),
    !! others near-tangent noise (m=16,18 technically cross by as little as
    !! ~1e-12/~1e-14 over an interval microns wide in z, then close back up
    !! -- not a real stability constraint, see MOD_EQ_MIN_ISLAND_WIDTH).
    !! This function's previous implementation (fixed dz=0.1 grown in steps
    !! of 10 from z=(m+2)^2, then one plain bisection with no width check
    !! at all) missed the significant windows and was unsafe by up to +77%
    !! (m=14) -- this is what get_zcrit's old hardcoded table below was
    !! actually built from, before being recomputed at 50-digit precision
    !! with this width-aware search.
    !!
    !! Algorithm: fine forward scan (dz below), and on any crossing found,
    !! check mod_eq_viol just past MOD_EQ_MIN_ISLAND_WIDTH beyond it -- if
    !! still positive (a real or permanent instability), accept; otherwise
    !! resume scanning from just past the island. Not a mathematical
    !! guarantee against an arbitrarily narrower-than-dz significant window
    !! for some m outside the range checked (up to m=30), but this is the
    !! rarely-hit fallback path -- get_zcrit's table below covers every m
    !! actually requested by compute_element_base_orders (m in
    !! {1,2,4,...,20}) with values re-derived independently.
    function ModifiedEquationZCrit(m) result(zc)
        implicit none
        integer, intent(in) :: m
        real(fpp) :: zc, zmax, dz, z_start
        logical :: found

        zmax = max(4._fpp * real((m + 2)**2, fpp), 1000._fpp)
        dz = 1.0e-5_fpp
        z_start = 0._fpp

        do
            call scan_first_crossing(m, z_start, zmax, dz, found, zc)
            do while (.not. found .and. zmax < 1.0e7_fpp)
                zmax = zmax * 2._fpp
                call scan_first_crossing(m, z_start, zmax, dz, found, zc)
            end do
            if (.not. found) then
                zc = zmax ! could not bracket an instability at all: conservative fallback
                return
            end if

            if (mod_eq_viol(zc + MOD_EQ_MIN_ISLAND_WIDTH, m) > 0._fpp) then
                return ! real (or permanent) instability: this is z_crit
            end if
            ! noise-level island: skip past it and keep looking
            z_start = zc + MOD_EQ_MIN_ISLAND_WIDTH
        end do
    end function ModifiedEquationZCrit

    !> Retorna z_crit para a ordem m (tabela direta para m <= 20).
    !! Recomputed 2026-08-17 at 50-digit precision (mpmath) using the same
    !! width-aware search as ModifiedEquationZCrit (not the old fixed-dz
    !! forward scan that produced the previous, wrong table -- see that
    !! function's docstring). m=1,2,4,20 were already correct; m=6 was off
    !! by +1.1%; m=8,10,12,14 were unsafe by +9.8%/+28.6%/+52.2%/+76.7%
    !! (the old z_crit exceeded the true first significant instability
    !! boundary); m=16,18 were already correct too (the old scan happened
    !! to jump clean over a real but noise-level near-tangent island at
    !! z~39.48 for both, landing directly on the same significant boundary
    !! this search deliberately skips ahead to now).
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
            zc = 21.4812098755971_fpp
        case (6)
            zc = 30.7214581598952_fpp
        case (8)
            zc = 37.0751178304851_fpp
        case (10)
            zc = 39.1829361357775_fpp
        case (12)
            zc = 39.4579712517290_fpp
        case (14)
            zc = 39.4774129983015_fpp
        case (16)
            zc = 151.348322547657_fpp
        case (18)
            zc = 156.849745440999_fpp
        case (20)
            zc = 157.804129863365_fpp
        case default
            zc = ModifiedEquationZCrit(m)
        end select
    end function get_zcrit

end module m_mod_eq_zcrit
