! Unit test for SolvePolynomial in mod_polynomialroots.F90, the quartic solver
! the discrete collision routines in mod_collisions.F90 depend on. Covers the
! bugs fixed in the July 2026 collision audit: the fourth root must be returned,
! the root pattern must be classified for the cubic and zero-constant paths,
! unused roots must come back as inert zeros, and outputCode must be
! threadprivate so concurrent solves from OpenMP regions do not corrupt it
! (the final case alternates two polynomials across threads to exercise this).
!
! Build against the compiled module object, e.g. from the build/ directory:
!   nvfortran -mp -O2 -I. ../Util/test_poly.f90 mod_polynomialroots.o -o test_poly.out
!   OMP_NUM_THREADS=16 ./test_poly.out
program test_poly
  use mod_polynomialroots
  implicit none
  complex(8) :: r1, r2, r3, r4
  integer :: code, i, nfail
  double precision :: re(4)
  integer :: nfail_omp

  nfail = 0

  ! Case 1: four real roots 1,2,3,4 -> code 31
  call SolvePolynomial(1.0d0, -10.0d0, 35.0d0, -50.0d0, 24.0d0, code, r1, r2, r3, r4)
  re = (/r1%re, r2%re, r3%re, r4%re/)
  call sort4(re)
  if (code /= 31) then; print *, 'FAIL case1 code', code; nfail = nfail + 1; end if
  if (any(abs(re - (/1.0d0,2.0d0,3.0d0,4.0d0/)) > 1.0d-8) .or. &
      abs(r1%im)+abs(r2%im)+abs(r3%im)+abs(r4%im) /= 0.0d0) then
    print *, 'FAIL case1 roots', re; nfail = nfail + 1
  end if

  ! Case 2: roots 1,2 and +/-i -> code 42, real roots first
  call SolvePolynomial(1.0d0, -3.0d0, 3.0d0, -3.0d0, 2.0d0, code, r1, r2, r3, r4)
  if (code /= 42) then; print *, 'FAIL case2 code', code; nfail = nfail + 1; end if
  if (r1%im /= 0.0d0 .or. r2%im /= 0.0d0) then; print *, 'FAIL case2 realroots'; nfail = nfail + 1; end if
  if (abs(min(r1%re,r2%re)-1.0d0) > 1d-8 .or. abs(max(r1%re,r2%re)-2.0d0) > 1d-8) then
    print *, 'FAIL case2 values', r1, r2; nfail = nfail + 1
  end if
  if (r3%im == 0.0d0 .or. r4%im == 0.0d0) then; print *, 'FAIL case2 compl'; nfail = nfail + 1; end if

  ! Case 3: all complex (t^2+1)(t^2+4) -> no root passes the "real and positive" scan
  call SolvePolynomial(1.0d0, 0.0d0, 5.0d0, 0.0d0, 4.0d0, code, r1, r2, r3, r4)
  if ((r1%im == 0 .and. r1%re > 0) .or. (r2%im == 0 .and. r2%re > 0) .or. &
      (r3%im == 0 .and. r3%re > 0) .or. (r4%im == 0 .and. r4%re > 0)) then
    if (code /= 44 .and. code /= 23) then
      print *, 'FAIL case3: complex root would pass scan, code', code; nfail = nfail + 1
    end if
  end if

  ! Case 4: zero constant term t(t-1)(t-2)(t-3) -> code must be set (31), roots 0,1,2,3
  call SolvePolynomial(1.0d0, -6.0d0, 11.0d0, -6.0d0, 0.0d0, code, r1, r2, r3, r4)
  re = (/r1%re, r2%re, r3%re, r4%re/)
  call sort4(re)
  if (code /= 31) then; print *, 'FAIL case4 code', code; nfail = nfail + 1; end if
  if (any(abs(re - (/0.0d0,1.0d0,2.0d0,3.0d0/)) > 1.0d-8)) then
    print *, 'FAIL case4 roots', re; nfail = nfail + 1
  end if

  ! Case 5: degenerate quadratic (acc=0): t^2-3t+2 -> code 22, root3/root4 inert zeros
  call SolvePolynomial(0.0d0, 0.0d0, 1.0d0, -3.0d0, 2.0d0, code, r1, r2, r3, r4)
  if (code /= 22) then; print *, 'FAIL case5 code', code; nfail = nfail + 1; end if
  if (r3 /= (0.0d0,0.0d0) .or. r4 /= (0.0d0,0.0d0)) then
    print *, 'FAIL case5 unused roots not zero', r3, r4; nfail = nfail + 1
  end if

  ! Case 6: threadprivate outputCode — alternate two polynomials across threads
  nfail_omp = 0
  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, code, r1, r2, r3, r4) REDUCTION(+:nfail_omp)
  do i = 1, 200000
    if (mod(i, 2) == 0) then
      call SolvePolynomial(1.0d0, -10.0d0, 35.0d0, -50.0d0, 24.0d0, code, r1, r2, r3, r4)
      if (code /= 31) nfail_omp = nfail_omp + 1
    else
      call SolvePolynomial(1.0d0, 0.0d0, 5.0d0, 0.0d0, 4.0d0, code, r1, r2, r3, r4)
      if (code == 31 .or. code == 42) nfail_omp = nfail_omp + 1
    end if
  end do
  !$OMP END PARALLEL DO
  if (nfail_omp /= 0) then; print *, 'FAIL case6 (race)', nfail_omp; nfail = nfail + 1; end if

  if (nfail == 0) then
    print *, 'ALL POLY TESTS PASSED'
  else
    print *, nfail, ' FAILURES'
  end if

contains
  subroutine sort4(x)
    double precision, intent(inout) :: x(4)
    integer :: p, q
    double precision :: tmp
    do p = 1, 3
      do q = p+1, 4
        if (x(q) < x(p)) then
          tmp = x(p); x(p) = x(q); x(q) = tmp
        end if
      end do
    end do
  end subroutine sort4
end program test_poly
