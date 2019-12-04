program elctrolyte

use mfkfun
use mkinsol

implicit none
integer i, j, counter,counter_max
real*8 z, ztrash
real*8, allocatable :: x_init(:)
character*22, densityfilename
character*13, mupolfilename
character*14, sysfilename

call readinput

neq = ntot 


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

Npol=Npol_min
counter=0
counter_max = int( (Npol_max - Npol_min)/Npol_step )

open(unit=4000,file='mupolvsnpol.dat')

do while (counter.le.counter_max)

  iter=0
  counter=counter+1
  
  write(sysfilename,'(A7,BZ,I3.3,A4)')'system.',counter,'.dat'
  write(densityfilename,'(A15,BZ,I3.3,A4)')'densitysolvent.',counter,'.dat'
  write(mupolfilename,'(A6,BZ,I3.3,A4)')'mupol.',counter,'.dat'

  open(unit=1000,file=densityfilename)
  open(unit=2000,file=mupolfilename)
  open(unit=3000,file=sysfilename)

  print*, "Solving Npol = ",Npol

  call call_kinsol(x_init)

  do i = 1, ntot
    z=(i-0.5)*delta
    write(1000,*)z, volumefraction(i)
    write(2000,*)z, mupol(i)

    x_init(i) = -log(volumefraction(i))

  enddo
  write(3000,*)"Npol= ", Npol
  write(3000,*)"st= ", st
  write(3000,*)"ntot= ", ntot
  write(3000,*)"vpol= ", vpol
  write(3000,*)"delta= ", delta
  write(3000,*)"lseg= ", lseg
 
  close(1000)
  close(2000)
  close(3000)
  
  write(4000,*)Npol,mupol(1)

  Npol=Npol+Npol_step
enddo

close(101)
close(4000)
end
