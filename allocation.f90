subroutine allocation

use mfkfun
use mkinsol

implicit none

allocate (osmoticpressure (ntot), volumefraction(ntot), mupol(ntot))
allocate (Xu(ntot,ntot))
allocate (pp(ntot))

end

