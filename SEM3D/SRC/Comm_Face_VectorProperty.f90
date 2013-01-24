subroutine Comm_Face_VectorProperty(ngll1,ngll2,orient_f,exch_val,prop_face)
    ! from vector of exchanged interproc values, to faces - vectorial case
    implicit none

    integer, intent(in)  :: ngll1,ngll2,orient_f
    real, intent(in), dimension(1:(ngll1-2)*(ngll2-2),0:2)  :: exch_val
    real, dimension(1:ngll1-2,1:ngll2-2,0:2), intent(inout) :: prop_face
    integer  :: j,k,ngll

    ngll = 1
    select case(orient_f)
    case(0)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(k,j,0:2) = prop_face(k,j,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(1)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-k,j,0:2) = prop_face(ngll1-1-k,j,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(2)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(k,ngll2-1-j,0:2) = prop_face(k,ngll2-1-j,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(3)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-k,ngll2-1-j,0:2) = prop_face(ngll1-1-k,ngll2-1-j,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(4)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(j,k,0:2) = prop_face(j,k,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(5)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-j,k,0:2) = prop_face(ngll1-1-j,k,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(6)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(j,ngll2-1-k,0:2) = prop_face(j,ngll2-1-k,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(7)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-j,ngll2-1-k,0:2) = prop_face(ngll1-1-j,ngll2-1-k,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    end select

    if(ngll /= (ngll1-2)*(ngll2-2)+1) &
        stop "Pb in Comm_Face_VectorProperty"

end subroutine Comm_Face_VectorProperty
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
