! Standalone timing harness: call lfmm3d_ndiv directly so we get the returned
! 6-phase timeinfo (1=P2M, 2=M2M, 3=M2L, 4=L2L, 5=L2P-eval, 6=P2P-direct) and
! tree+setup = total - sum(timeinfo). Shows how each phase scales with nd, and
! that tree construction is a tiny, nd-independent fraction.
!
! Build (needs an FMM3D checkout with the static lib already built; set FMM3D):
!   FMM3D=$HOME/codes/FMM3D
!   gfortran -O3 -fopenmp -march=native nd_timing.f90 \
!       $FMM3D/lib-static/libfmm3d.a -o nd_timing -lgomp -lm -lstdc++ -lgfortran -ldl
! Run:  OMP_NUM_THREADS=1 ./nd_timing
!
! Result (N=200k uniform, 1 thread, eps=1e-4): tree+setup ~2% of the FMM; the
! shared ~84% is per-interaction geometry/operators (M2L + P2P + Legendre/1/r)
! reused across all nd via the `do idim=1,nd` inner loop in every kernel.
program nd_timing
  implicit none
  integer*8 :: ns, nt, nd, ifcharge, ifdipole, iper, ifpgh
  integer*8 :: ntarg, ifpghtarg, ifnear, ndiv, idivflag, ier
  integer*8 :: i, idim, k
  integer*8 :: ndlist(2)
  real*8 :: eps, t0, t1, ttot, tph, omp_get_wtime, timeinfo(6)
  real*8 :: dipvec(3,1), targ(3,1), pottarg(1,1), gradtarg(1,3,1), hesstarg(1,6,1)
  real*8, allocatable :: source(:,:), charge(:,:), pot(:,:), grad(:,:,:), hess(:,:,:)
  external omp_get_wtime

  ns = 200000; nt = 0; eps = 1.0d-4
  ifcharge = 1; ifdipole = 0; iper = 0; ifpgh = 2
  ntarg = 0; ifpghtarg = 0; ifnear = 1
  ndlist = (/ 1_8, 36_8 /)

  allocate(source(3,ns))
  call random_seed()
  do i = 1, ns
     call random_number(source(1,i))
     call random_number(source(2,i))
     call random_number(source(3,i))
  end do

  do k = 1, 2
     nd = ndlist(k)
     allocate(charge(nd,ns), pot(nd,ns), grad(nd,3,ns), hess(nd,6,ns))
     do i = 1, ns
        do idim = 1, nd
           call random_number(charge(idim,i))
        end do
     end do
     pot = 0; grad = 0; hess = 0
     call lndiv(eps, ns, nt, ifcharge, ifdipole, ifpgh, ifpghtarg, ndiv, idivflag)
     t0 = omp_get_wtime()
     call lfmm3d_ndiv(nd, eps, ns, source, ifcharge, charge, ifdipole, dipvec, &
          iper, ifpgh, pot, grad, hess, ntarg, targ, ifpghtarg, pottarg, &
          gradtarg, hesstarg, ndiv, idivflag, ifnear, timeinfo, ier)
     t1 = omp_get_wtime()
     ttot = t1 - t0
     tph  = sum(timeinfo)
     write(*,'(A,I3,A,I6,A,I2)') '==== nd=', nd, '  ndiv=', ndiv, '  ier=', ier
     write(*,'(A,F9.3)') '  total          ', ttot
     write(*,'(A,F9.3)') '  P2M  step1     ', timeinfo(1)
     write(*,'(A,F9.3)') '  M2M  step2     ', timeinfo(2)
     write(*,'(A,F9.3)') '  M2L  step3     ', timeinfo(3)
     write(*,'(A,F9.3)') '  L2L  step4     ', timeinfo(4)
     write(*,'(A,F9.3)') '  L2P  step5     ', timeinfo(5)
     write(*,'(A,F9.3)') '  P2P  step6     ', timeinfo(6)
     write(*,'(A,F9.3)') '  sum(phases)    ', tph
     write(*,'(A,F9.3)') '  tree+setup     ', ttot - tph
     deallocate(charge, pot, grad, hess)
  end do
end program nd_timing
