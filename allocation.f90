subroutine allocation

use mfkfun
use mkinsol

implicit none

allocate (osmoticpressure (ntot), volumefraction(ntot))
allocate (Xu(ntot,ntot))
allocate (pp(neq))

end

