module mkinsol

integer*8 neq
real*8, allocatable :: pp(:)

endmodule mkinsol

module mfkfun

real*8 Npol,vpol,st,lseg !caracteriticas de los bichos
integer*8 ntot, Xulimit, flagkai
real*8, allocatable :: osmoticpressure(:), volumefraction(:), Xu(:,:)
real*8, parameter :: NA=6.02d23, Eps=0.114, Beta=1 !Constantes. Eps tiene unidades de e^2/kT.nm
real*8 delta  !layer length (nm)


endmodule mfkfun
