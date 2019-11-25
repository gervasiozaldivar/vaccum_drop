program elctrolyte

use mfkfun
use mkinsol

implicit none
integer i, j
real*8 z, ztrash
real*8, allocatable :: x_init(:)

call readinput

neq = ntot*2 
iter=0


call allocation


print*,"kai calculation"
call kai

allocate (x_init(ntot*2))

x_init = 1.0 ! homogeneous initial guess

if (infile.eq.1) then
  do i=1,ntot
    read(101,*), ztrash, x_init(i)
    read(102,*), ztrash, x_init(ntot+i)
  enddo
  infile = 2
endif

call call_kinsol(x_init)

do i = 1, ntot
 z=i*delta
 write(101,*)z, volumefraction(i)
enddo

close(101)

end
