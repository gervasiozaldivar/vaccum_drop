subroutine fkfun (x, f, ier) !x, f y ier son entradas y/o salidas de la subrutina, para comunicarse con otras subrutinas

use mfkfun

implicit none

integer*8 j,i !indices
integer ier
real*8 f(*), x(*)  ! el asterisco es porque x y f entran desde y salen al kinsol. Es decirle que use el tamaño que está indicado en kinsol para f. 
                   ! x es rho(s)*vol(s) 

Cbulk(1)=(1d24/NA-Cbulk(2)*vol(2) - Cbulk(3)*vol(3))/vol(1)

psi(1:ntot) = x(ntot+1:2*ntot)
psi(0) = psi(1)+sigma*delta/Eps !condiciones de contorno para Psi(0)
psi(ntot+1) = 0.0 !condiciones de contorno para Psi(bulk)

f(1:2*ntot) = 0.0

print*, 'fkfun entra al loop'

do i=1,ntot

     do j=1,3

     rho(i,j)=Cbulk(j)*NA*1e-24/( (Cbulk(1)*vol(1)*NA*1e-24)**(vol(j)/vol(1)))*(x(i)**(vol(j)/vol(1)))*exp(Beta*(-Psi(i)*z(j))) !está escrito para x(1:ntot) = rho(s)*vol(s)
!     rho(i,j)=Cbulk(j)*NA*1e-24*exp(Beta*(-pi(i)*vol(j)/vol(1)-Psi(i)*z(j))) ! está escrito para x(1:ntot) = pi

     f(i) = f(i) + vol(j)*rho(i,j)

     f(ntot + i) = f(ntot + i) + z(j)*rho(i,j)/Eps

     enddo

f(i) = f(i) - 1
f(ntot + i) = f(ntot + i) + (Psi(i+1)-2*Psi(i)+Psi(i-1))*delta**(-2)
f(ntot + i) = f(ntot + i)/2.0 !por qué esto?

enddo

print*, 'fkfun deja de iterar'

ier = 0 !si ier ne 0 el kinsol tira error.

end

