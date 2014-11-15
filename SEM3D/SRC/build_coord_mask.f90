subroutine build_coord_mask(Tdomain, mat)
    use sdomain

    implicit none

    type(domain), intent(inout), target :: Tdomain
    integer              :: ngllx, nglly, ngllz
    integer              :: m, n, i, j, k, coord, mat, ipoint

	if(.not.allocated(Tdomain%sSubDomain(mat)%globCoordMask))                          &
		allocate(Tdomain%sSubDomain(mat)%globCoordMask(0:size(Tdomain%GlobCoord,1)-1,  &
	        	 									   0:size(Tdomain%GlobCoord,2)-1))
    Tdomain%sSubDomain(mat)%globCoordMask(:,:) = .false.
    !write(*,*)"...%globCoordMask(:,:) = BEFORE", Tdomain%sSubDomain(mat)%globCoordMask(:,:)

    ngllx = Tdomain%sSubDomain(mat)%NGLLx
    nglly = Tdomain%sSubDomain(mat)%NGLLy
    ngllz = Tdomain%sSubDomain(mat)%NGLLz

    do m = 0, Tdomain%sSubDomain(mat)%nElem - 1
        n = Tdomain%sSubDomain(mat)%elemList(m)
        do i = 0, ngllx-1
            do j = 0, nglly-1
                do k = 0, ngllz-1
                    ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)
                    Tdomain%sSubDomain(mat)%globCoordMask(:,ipoint) = .true.
                end do
            end do
    	end do !END Loop over GLLs
    end do !END Loop over subdomain elements
	!write(*,*)"...%globCoordMask(:,:) AFTER = ", Tdomain%sSubDomain(mat)%globCoordMask(:,:)
end subroutine build_coord_mask
