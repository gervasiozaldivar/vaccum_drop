program elctrolyte

use mfkfun
use mkinsol

implicit none
integer i, j
real*8 z
call readinput

neq = ntot*2 
iter=0


call allocation


print*,"kai calculation"
call kai

call call_kinsol

do i = 1, ntot
 z=i*delta
 write(101,*)z, volumefraction(i)
enddo

close(101)

end
