subroutine fkfun (x, f, ier) !x, f y ier son entradas y/o salidas de la subrutina, para comunicarse con otras subrutinas

use mfkfun

implicit none

integer*8 j,i !indices
integer ier
real*8 f(ntot), x(ntot)      ! x(1:ntot)=volumefraction(i) / x(ntot+1,ntot*2)=pi(i) 
real*8 suminteractions, repulsions, sumpol, packing, exponente
real*8 pi_kinsol(ntot),volumefraction_kinsol (ntot)
real*8 algo

iter=iter+1

! pi_kinsol(1:ntot) = x(ntot+1:2*ntot) ! pi is read from kinsol x
volumefraction_kinsol(1:ntot) = x(1:ntot) ! volume fraction is read from kinsol x

f(1:ntot) = 0.0


sumpol=0.0

do i=1,ntot

suminteractions=0.0

  do j=1,ntot
    suminteractions = suminteractions + Xu(i,j)*st*volumefraction_kinsol(j) 
  enddo
  
  repulsions = 8.0*volumefraction_kinsol(i) - 9.0*volumefraction_kinsol(i) + 3.0*volumefraction_kinsol(i) 
  repulsions = repulsions * ( 1 - volumefraction_kinsol (i) ) ** (-3.0)

  volumefraction(i) =  exp(suminteractions-repulsions)

  sumpol = sumpol + volumefraction(i)*delta ! 
  

enddo

volumefraction = volumefraction * Npol * vpol / sumpol   


do i=1,ntot
  f(i) = volumefraction(i)-volumefraction_kinsol(i)
enddo


ier = 0 !si ier ne 0 el kinsol tira error.

algo=0.0
do i=1,ntot*2
  algo=algo+f(i)**2
enddo

print*,iter,algo

end

