program elctrolyte

use mfkfun
use mkinsol

implicit none
integer i, j
real*8 z, ztrash
real*8, allocatable :: x_init(:)

call readinput

neq = ntot 
iter=0


call allocation


print*,"kai calculation"

call kai

allocate (x_init(ntot))

x_init = 0.1 ! homogeneous initial guess

if (infile.eq.1) then
  do i=1,ntot
    read(101,*), ztrash, x_init(i)
  enddo
  infile = 2
endif

call call_kinsol(x_init)

do i = 1, ntot
  z=(i-0.5)*delta
  write(1000,*)z, volumefraction(i)
  write(2000,*)z, mupol(i)
enddo

close(101)

end
