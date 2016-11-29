!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module tensor_util
    use constants, only : fpp
    implicit none
    real(fpp), parameter :: s2=1.4142135623730950488d0
    real(fpp), parameter :: os2=0.70710678118654d0

contains
    ! #########################################################
    subroutine c_4tensor(CIJ,theta,phi)
    ! transformation de spherique a cartesien
        implicit none

        real(fpp), dimension(1:6,1:6), intent(inout) :: CIJ
        real(fpp), intent(in) :: theta,phi

        integer :: i,j,k,l,a,b,c,d
        real(fpp) :: ct,st,cp,sp
        real(fpp), dimension(3,3) :: CCS
        real(fpp), dimension(3,3,3,3) :: cijkl,T

        call kelvin2ijkl(CIJ,T)

        ! passage de cartersien a spherique
        ct=dcos(theta)
        st=dsin(theta)
        cp=dcos(phi  )
        sp=dsin(phi  )
        Ccs(1,1)= st*cp; Ccs(1,2)=st*sp; Ccs(1,3)=  ct
        Ccs(2,1)= ct*cp; Ccs(2,2)=ct*sp; Ccs(2,3)= -st
        Ccs(3,1)=-sp   ; Ccs(3,2)=cp   ; Ccs(3,3)= 0.d0

        ! M*Ccs
        do i=1,3
            do j=1,3
                do k=1,3
                    do l=1,3
                        cijkl(i,j,k,l)=0.d0
                        do a=1,3
                            do b=1,3
                                do c=1,3
                                    do d=1,3
                                        cijkl(i,j,k,l)=cijkl(i,j,k,l)+T(a,b,c,d)*Ccs(a,i)*Ccs(b,j)*Ccs(c,k)*Ccs(d,l)
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        call ijkl2kelvin(cijkl,CIJ)

    end subroutine c_4tensor

    ! #########################################################
    subroutine rot_4tensor(CIJ,rot)

        ! rotation d'un chunk a l'autre

        implicit none

        real(fpp), dimension(1:6,1:6), intent(inout) :: CIJ
        real(fpp), dimension(1:3,1:3), intent(in) :: rot

        integer :: i,j,k,l,a,b,c,d
        real(fpp), dimension(3,3,3,3) :: cijkl,T

        call kelvin2ijkl(CIJ,T)

        do i=1,3
            do j=1,3
                do k=1,3
                    do l=1,3    
                        cijkl(i,j,k,l)=0.d0
                        do a=1,3
                            do b=1,3
                                do c=1,3
                                    do d=1,3
                                        cijkl(i,j,k,l)=cijkl(i,j,k,l)+T(a,b,c,d)*rot(a,i)*rot(b,j)*rot(c,k)*rot(d,l)
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo 
                enddo
            enddo 
        enddo
        call ijkl2kelvin(cijkl,CIJ)

        end subroutine rot_4tensor

! #########################################################
    subroutine ijkl2kelvin(cijkl,CIJ)

        implicit none

        real(fpp), dimension(1:6,1:6), intent(out) :: CIJ
        real(fpp), dimension(1:3,1:3,1:3,1:3), intent(in) :: cijkl

        integer :: i,j

        CIJ(1,1)=     cijkl(1,1,1,1)
        CIJ(2,2)=     cijkl(2,2,2,2)
        CIJ(3,3)=     cijkl(3,3,3,3)
        CIJ(1,2)=     cijkl(1,1,2,2)
        CIJ(1,3)=     cijkl(1,1,3,3)
        CIJ(2,3)=     cijkl(2,2,3,3)

        CIJ(4,4)=2.d0*cijkl(2,3,2,3)
        CIJ(5,5)=2.d0*cijkl(1,3,1,3)
        CIJ(6,6)=2.d0*cijkl(1,2,1,2)
        CIJ(4,5)=2.d0*cijkl(2,3,1,3)
        CIJ(4,6)=2.d0*cijkl(2,3,1,2)
        CIJ(5,6)=2.d0*cijkl(1,3,1,2)

        CIJ(1,4)=  s2*cijkl(1,1,2,3)
        CIJ(1,5)=  s2*cijkl(1,1,1,3)
        CIJ(1,6)=  s2*cijkl(1,1,1,2)
        CIJ(2,4)=  s2*cijkl(2,2,2,3)
        CIJ(2,5)=  s2*cijkl(2,2,1,3)
        CIJ(2,6)=  s2*cijkl(2,2,1,2)
        CIJ(3,4)=  s2*cijkl(3,3,2,3)
        CIJ(3,5)=  s2*cijkl(3,3,1,3)
        CIJ(3,6)=  s2*cijkl(3,3,1,2)

        do i=2,6 
            do j=1,i-1
                CIJ(i,j)=CIJ(j,i)
            enddo
        enddo

    end subroutine ijkl2kelvin

    ! ##########################################################
    subroutine kelvin2ijkl(CIJ,cijkl)

        implicit none

        real(fpp), dimension(1:6,1:6), intent(in) :: CIJ
        real(fpp), dimension(1:3,1:3,1:3,1:3), intent(out) :: cijkl

        integer :: i,j,k,l,ii
        integer, dimension(4,21) :: ijkltab

        call init_ijkltab(ijkltab)

        cijkl(1,1,1,1)=      CIJ(1,1)
        cijkl(2,2,2,2)=      CIJ(2,2)
        cijkl(3,3,3,3)=      CIJ(3,3)
        cijkl(1,1,2,2)=      CIJ(1,2)
        cijkl(1,1,3,3)=      CIJ(1,3)
        cijkl(2,2,3,3)=      CIJ(2,3)

        cijkl(2,3,2,3)=0.5d0*CIJ(4,4)
        cijkl(1,3,1,3)=0.5d0*CIJ(5,5)
        cijkl(1,2,1,2)=0.5d0*CIJ(6,6)
        cijkl(2,3,1,3)=0.5d0*CIJ(4,5)
        cijkl(2,3,1,2)=0.5d0*CIJ(4,6)
        cijkl(1,3,1,2)=0.5d0*CIJ(5,6)

        cijkl(1,1,2,3)=  os2*CIJ(1,4)
        cijkl(1,1,1,3)=  os2*CIJ(1,5)
        cijkl(1,1,1,2)=  os2*CIJ(1,6)
        cijkl(2,2,2,3)=  os2*CIJ(2,4)
        cijkl(2,2,1,3)=  os2*CIJ(2,5)
        cijkl(2,2,1,2)=  os2*CIJ(2,6)
        cijkl(3,3,2,3)=  os2*CIJ(3,4)
        cijkl(3,3,1,3)=  os2*CIJ(3,5)
        cijkl(3,3,1,2)=  os2*CIJ(3,6)

       do ii=1,21
          i=ijkltab(1,ii)
          j=ijkltab(2,ii)
          k=ijkltab(3,ii)
          l=ijkltab(4,ii)

          cijkl(j,i,k,l)= cijkl(i,j,k,l)
          cijkl(i,j,l,k)= cijkl(i,j,k,l)
          cijkl(j,i,l,k)= cijkl(i,j,k,l)

          cijkl(k,l,i,j)= cijkl(i,j,k,l)
          cijkl(l,k,i,j)= cijkl(i,j,k,l)
          cijkl(k,l,j,i)= cijkl(i,j,k,l)
          cijkl(l,k,j,i)= cijkl(i,j,k,l)   
       enddo

       contains
!---------------------------------------
       subroutine init_ijkltab(ijkltab)

         implicit none
         integer, dimension(:,:), intent(out) :: ijkltab

         ijkltab(:,1 )=(/1,1,1,1/)
         ijkltab(:,2 )=(/2,2,2,2/)
         ijkltab(:,3 )=(/3,3,3,3/)
         ijkltab(:,4 )=(/1,1,2,2/)
         ijkltab(:,5 )=(/1,1,3,3/)
         ijkltab(:,6 )=(/2,2,3,3/)

         ijkltab(:,7 )=(/2,3,2,3/)
         ijkltab(:,8 )=(/1,3,1,3/)
         ijkltab(:,9 )=(/1,2,1,2/)
         ijkltab(:,10)=(/2,3,1,3/)
         ijkltab(:,11)=(/2,3,1,2/)
         ijkltab(:,12)=(/1,3,1,2/)

         ijkltab(:,13)=(/1,1,2,3/)
         ijkltab(:,14)=(/1,1,1,3/)
         ijkltab(:,15)=(/1,1,1,2/)
         ijkltab(:,16)=(/2,2,2,3/)
         ijkltab(:,17)=(/2,2,1,3/)
         ijkltab(:,18)=(/2,2,1,2/)
         ijkltab(:,19)=(/3,3,2,3/)
         ijkltab(:,20)=(/3,3,1,3/)
         ijkltab(:,21)=(/3,3,1,2/)

       end subroutine init_ijkltab
!---------------------------------------

    end subroutine kelvin2ijkl

! ###############################################
    real(fpp) function lambda_from_Cij (C)

        real(fpp), dimension(1:6,1:6), intent(IN) :: C

        lambda_from_Cij = ( C(1,1) + C(2,2) + C(3,3) &
                            + 4.d0*C(1,2) + 4.d0*C(2,3) + 4.d0*C(1,3) &
                            - C(4,4) - C(5,5) - C(6,6) ) / 15.d0

        return
    end function

    ! ###############################################
    real(fpp) function mu_from_Cij (C)

        real(fpp), dimension(1:6,1:6), intent(IN) :: C

        mu_from_Cij = ( C(1,1) + C(2,2) + C(3,3) &
                      - C(1,2) - C(2,3) - C(1,3) &
                      + 1.5*C(4,4) + 1.5*C(5,5) + 1.5*C(6,6) ) / 15.d0

        return
    end function

! ############################################################
end module tensor_util

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
