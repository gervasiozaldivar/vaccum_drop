subroutine fkfun (x, f, ier) !x, f y ier son entradas y/o salidas de la subrutina, para comunicarse con otras subrutinas

use mfkfun

implicit none

integer*8 j,i !indices
integer ier
real*8 f(ntot), x(ntot)      ! x(1:ntot)=volumefraction(i)  
real*16 suminteractions(ntot), repulsions(ntot), sumpol, packing, exponente
real*16 pi_kinsol(ntot),volumefraction_kinsol(ntot)
real*16 algo, denominador, st_pared(ntot)

iter=iter+1

! pi_kinsol(1:ntot) = x(ntot+1:2*ntot) ! pi is read from kinsol x
do i =1,ntot
  volumefraction_kinsol(i) = exp(-x(i)) ! volume fraction is read from kinsol x
enddo

f(1:ntot) = 0.0

sumpol=0.0
st_pared=0.0
!st_pared(1)=-5.0


  suminteractions = 0.0
  repulsions = 0.0

do i=1,ntot

  denominador = 0.0

    do j=1,ntot
      suminteractions(i) = suminteractions(i) + Xu(i,j)*st*volumefraction_kinsol(j)/vpol 
    enddo
  
  denominador = (1.0 - volumefraction_kinsol(i))
  repulsions(i)=(8.0*volumefraction_kinsol(i)-9.0*volumefraction_kinsol(i)**2+3.0*volumefraction_kinsol(i)**3) / denominador**3 

  volumefraction(i) =  exp(suminteractions(i)-repulsions(i)-st_pared(i))
  
  
  sumpol = sumpol + volumefraction(i)*delta ! 

enddo

volumefraction = volumefraction * Npol * vpol / sumpol   

do i=1,ntot
  f(i) = volumefraction(i)-volumefraction_kinsol(i)
  mupol(i) = log(volumefraction(i)) + repulsions(i) - suminteractions(i) - st_pared(i)
enddo


ier = 0 !si ier ne 0 el kinsol tira error.

algo=0.0
do i=1,ntot
  algo=algo+f(i)**2
!  print*, i, f(i), x(i)
enddo

print*,iter,algo

end

