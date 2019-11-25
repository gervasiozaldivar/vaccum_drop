program elctrolyte

use mfkfun
use mkinsol

implicit none
integer i, j

call readinput

neq = ntot*2 

call allocation

call call_kinsol

do i = 1, ntot
 write(100+j,*)i, rho(i,j)
enddo

close(101)

end
