subroutine fkfun (x, f, ier) !x, f y ier son entradas y/o salidas de la subrutina, para comunicarse con otras subrutinas

use mfkfun

implicit none

integer*8 j,i !indices
integer ier
real*8 f(ntot*2), x(ntot*2)      ! x(1:ntot)=volumefraction(i) / x(ntot+1,ntot*2)=pi(i) 
real*8 suminteractions, sumpol, packing, exponente
real*8 pi_kinsol(ntot),volumefraction_kinsol (ntot)
real*8 algo

iter=iter+1

pi_kinsol(1:ntot) = x(ntot+1:2*ntot) ! pi is read from kinsol x
volumefraction_kinsol(1:ntot) = x(1:ntot) ! volume fraction is read from kinsol x

f(1:2*ntot) = 0.0


sumpol=0.0

do i=1,ntot
if (iter.eq.1) write(200,*), i, volumefraction_kinsol(i)
suminteractions=0.0

  do j=1,ntot
    suminteractions = suminteractions + Xu(i,j)*st*volumefraction_kinsol(j) 
  enddo

  exponente = -pi_kinsol(i)*vpol + suminteractions

  volumefraction(i) =  exp(exponente)

  sumpol = sumpol + volumefraction(i)*delta ! 
  

enddo

volumefraction = volumefraction * Npol * vpol / sumpol   

do i=1,ntot
  f(i) = volumefraction(i)-volumefraction_kinsol(i)
  packing = volumefraction(i)-1.0
  if (packing.lt.0.0) then
    f(ntot+i) = pi_kinsol(i)
    else
      f(ntot+i) = packing
  endif 
enddo


ier = 0 !si ier ne 0 el kinsol tira error.

algo=0.0
do i=1,ntot*2
  algo=algo+f(i)**2
enddo

print*,iter,algo

end

