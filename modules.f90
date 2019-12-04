module mkinsol

integer*8 neq
real*8, allocatable :: pp(:)

endmodule mkinsol

module mfkfun

real*16 Npol,vpol,st,lseg !caracteriticas de los bichos
integer*8 ntot, Xulimit, flagkai, infile
real*16 Npol_min, Npol_max, Npol_step
real*16, allocatable :: osmoticpressure(:), volumefraction(:), Xu(:,:), mupol(:)
real*16, parameter :: NA=6.02d23, Eps=0.114, Beta=1 !Constantes. Eps tiene unidades de e^2/kT.nm
real*16 delta  !layer length (nm)
integer*8 iter

endmodule mfkfun
