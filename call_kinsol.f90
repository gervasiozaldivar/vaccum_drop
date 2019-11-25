subroutine call_kinsol(x_in)
      use mkinsol 
      use mfkfun
      implicit none
      
      integer*8 max_niters
      integer*8 iout(15), msbpre ! bichos de kinsol
      integer ier, globalstrat, maxl, maxlrst, i ! bichos de kinsol
      double precision fnormtol, scsteptol 
      double precision rout(2)
      real*8, allocatable :: xg(:), scale(:), constr(:)
      integer h
      real*8 x_in(ntot*2)

      !ntot = 10
      !neq = ntot

      allocate (xg(neq))
      allocate (scale(neq))
      allocate (constr(neq))
      !allocate (pp(neq))

      globalstrat = 0 !el kinsol tiene varias capas de aproximación, esta es la primera capa. 0/1 indica tipo de solver.
      fnormtol = 1.0d-8 !tolerancia de la norma de f (debe ser cero)
      scsteptol = 1.0d-8 !tolerancia del cambio de f entre iteración e iteración
      maxl = 2000 !Krislov
      maxlrst = 50 !?
      msbpre  = 10 !?
      max_niters = 1000 !max. numero de iteraciones
!c * * * * * * * * * * * * * * * * * * * * * *

      call fnvinits(3, neq, ier) !kinsol setea neq. come neq, vomita ier

      if (ier .ne. 0) then
         write(6,1220) ier
 1220    format('SUNDIALS_ERROR: FNVINITS returned IER = ', i4)
         stop
      endif

!Defino un initial guess tipo bulk

      do i = 1, neq  ! 
         xg(i) = x_in(i) ! read from main.f90 
         scale(i) = 1.0d0 ! le da mas o menos importancia a cada variable a la hora de calcular el residual
         constr(i) = 2.0d0 ! restricciones para x, puede tomar 5 valores, que indican positivo, negativo, positivo esctricto, negativo estricto, o cualquier cosa (0)
      enddo

 
      !xg(5) = 1.0

      call fkincreate(ier)
      if (ier .ne. 0) then
         write(6,1230) ier
 1230    format('SUNDIALS_ERROR: FKINCREATE returned IER = ', i4)
         stop
      endif

      call fkinsetiin('MAX_SETUPS', msbpre, ier)
      if (ier .ne. 0) then
         write(6,1231) ier
 1231    format('SUNDIALS_ERROR: FKINSETIIN returned IER = ', i4)
         call fkinfree
         stop
      endif
      
      call fkinsetiin('MAX_NITERS', max_niters, ier) !las mayus entre comillas son lo que va leer kinsol(C) y la primer variable es lo que le vamos con ese nombre
      if (ier .ne. 0) then
         write(6,1231) ier
 1241    format('SUNDIALS_ERROR: FKINSETIIN returned IER = ', i4)
         call fkinfree
         stop
      endif


      call fkinsetrin('FNORM_TOL', fnormtol, ier)
      if (ier .ne. 0) then
         write(6,1232) ier
 1232    format('SUNDIALS_ERROR: FKINSETRIN returned IER = ', i4)
         call fkinfree
         stop
      endif

      call fkinsetrin('SSTEP_TOL', scsteptol, ier)
      if (ier .ne. 0) then
         write(6,1232) ier
         call fkinfree
         stop
      endif

      call fkinsetvin('CONSTR_VEC', constr, ier)
      if (ier .ne. 0) then
         write(6,1233) ier
 1233    format('SUNDIALS_ERROR: FKINSETVIN returned IER = ', i4)
         call fkinfree
         stop
      endif

      call fkininit(iout, rout, ier)
      if (ier .ne. 0) then
         write(6,1234) ier
 1234    format('SUNDIALS_ERROR: FKININIT returned IER = ', i4)
         stop
      endif

      call fkinspgmr(maxl, maxlrst, ier)
      if (ier .ne. 0) then
         write(6,1235) ier
 1235    format('SUNDIALS_ERROR: FKINSPGMR returned IER = ', i4)
         call fkinfree
         stop
      endif

      call fkinspilssetprec(1, ier)

!      write(6,1240)
! 1240 format('Example program fkinDiagon_kry:'//' This FKINSOL example'&
!     1       ' solves a 128 eqn diagonal algebraic system.'/
!     2       ' Its purpose is to demonstrate the use of the Fortran'&
!     3       ' interface'/' in a serial environment.'///
!     4       ' globalstrategy = KIN_NONE')


      print*, 'llamamos a fkinsol'


!------------------------llamamos a fkinsol

      call fkinsol(xg, globalstrat, scale, scale, ier)
!-----------------------kinsol escupe x
      print*, 'kinsol escupe x'

!      do h=1,ntot
!      print*, rho(h,:)
!      enddo

      if (ier .lt. 0) then
         write(6,1242) ier, iout(9)
 1242    format('SUNDIALS_ERROR: FKINSOL returned IER = ', i4, '     Linear Solver returned IER = ', i4)
         call fkinfree
         stop
      endif


      write(6,1245) ier
 1245 format(/' FKINSOL return code is ', i4)

      write(6,1246)
 1246 format(//' The resultant values of xg are:'/)

      do 30 i = 1, neq, 4
         write(6,1256) i, xg(i), xg(i+1), xg(i+2), xg(i+3)
 1256    format(i4, 4(1x, f10.6))
 30   continue

!      write(6,1267) iout(3), iout(14), iout(4), iout(12), iout(13), iout(15)

! 1267 format(//'Final statistics:'//
!          '  nni = ', i3, ',  nli  = ', i3, /,
!          '  nfe = ', i3, ',  npe  = ', i3, /,
!          '  nps = ', i3, ',  ncfl = ', i3)


      call fkinfree

!     stop
      end
      
! ------------Despedimos a call_kinsol------------

!c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!c     The routine kpreco is the preconditioner setup routine. It must have
!c     that specific name be used in order that the c code can find and link
!c     to it.  The argument list must also be as illustrated below:
      
      subroutine fkpset(udata, uscale, fdata, fscale, vtemp1, vtemp2, ier)
      
!subroutine fkpset(udata, uscale, vtemp1, vtemp2, ier)

      use mkinsol

      implicit none

      integer ier, i
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*)

      !allocate (pp(neq))

      do 10 i = 1, neq
         pp(i) = 1.0 !/ (udata(i) + 5.0d0) 
 10   continue
      ier = 0
      
      return
      end
      
      
!c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!c     The routine kpsol is the preconditioner solve routine. It must have
!c     that specific name be used in order that the c code can find and link
!c     to it.  The argument list must also be as illustrated below:
      
      subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)
      !subroutine fkpsol(vv, ier)
      use mkinsol

      implicit none

      integer ier, i

      double precision udata(*), uscale(*), fdata(*), fscale(*)
!      double precision udata(*), uscale(*)
      double precision vv(*), ftem(*)
      !double precision vv(*)

      !common /pcom/ pp(128)
      !common /psize/ neq

      !allocate (pp(neq))

      do 10 i = 1, neq
         vv(i) = vv(i) * pp(i)
 10   continue
      ier = 0
      
      return
      end
