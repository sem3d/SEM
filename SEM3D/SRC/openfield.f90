subroutine openfield (Tdomain,rg)

! Modified by Elise Delavaud (14/03/06)

use sdomain

implicit none

type (domain), intent (IN):: Tdomain
integer, intent (IN) :: rg

! local variables
character (len=20) :: fnamef

write (fnamef,"(a,I2.2,a)") "Proc",rg,".flavia.res"
open (61, file=fnamef, status="unknown", form="formatted")
write (61,*) 'GiD Post Results File 1.0'
write (61,*) 
write (61,*) 'ResultRangesTable "My Table"'
write (61,*) '0.0 - 0.5: "blabla1"'
write (61,*) '0.5 - 1.0: "blabla2"'
write (61,*) 'End ResultRangesTable'
write (61,*) 


return
end subroutine openfield
