! Benchmark: pack-based vs manual in-place compaction, as used in
! mod_pair.F90 Remove_Particles. Mirrors the real call pattern:
! double 1D, integer 1D, and (3,N) double arrays, restricted to 1:m.
program bench_compact
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: m, ndead, rep, nrep, scenario
  double precision, allocatable :: a1(:), a1ref(:), a2(:,:), a2ref(:,:)
  integer, allocatable :: ai(:), airef(:)
  logical, allocatable :: mask(:)
  double precision :: t0, t1, tpack1, tman1, tpack2, tman2, tpacki, tmani
  integer :: k0, i

  do scenario = 1, 3
    select case (scenario)
    case (1) ! large system, 1% scattered removals
      m = 1000000; ndead = 10000; nrep = 50
      print '(a)', '--- 1,000,000 particles, 10,000 scattered removals ---'
    case (2) ! large system, few removals
      m = 1000000; ndead = 10; nrep = 50
      print '(a)', '--- 1,000,000 particles, 10 scattered removals ---'
    case (3) ! small system, few removals
      m = 10000; ndead = 20; nrep = 5000
      print '(a)', '--- 10,000 particles, 20 scattered removals ---'
    end select

    allocate(a1(m), a1ref(m), a2(3,m), a2ref(3,m), ai(m), airef(m), mask(m))
    call random_number(a1ref)
    call random_number(a2ref)
    airef = [(i, i=1,m)]
    call make_mask(mask, m, ndead)
    k0 = findloc(mask, .false., dim=1)

    ! 1D double
    tpack1 = 0.0d0; tman1 = 0.0d0
    do rep = 1, nrep
      a1 = a1ref
      call cpu_time(t0)
      a1(1:m) = pack(a1(1:m), mask(1:m), a1(1:m))
      call cpu_time(t1)
      tpack1 = tpack1 + (t1-t0)
    end do
    do rep = 1, nrep
      a1 = a1ref
      call cpu_time(t0)
      call man1(a1, mask, k0, m)
      call cpu_time(t1)
      tman1 = tman1 + (t1-t0)
    end do

    ! 2D double (3,N), component-wise pack as in mod_pair
    tpack2 = 0.0d0; tman2 = 0.0d0
    do rep = 1, nrep
      a2 = a2ref
      call cpu_time(t0)
      a2(1, 1:m) = pack(a2(1, 1:m), mask(1:m), a2(1, 1:m))
      a2(2, 1:m) = pack(a2(2, 1:m), mask(1:m), a2(2, 1:m))
      a2(3, 1:m) = pack(a2(3, 1:m), mask(1:m), a2(3, 1:m))
      call cpu_time(t1)
      tpack2 = tpack2 + (t1-t0)
    end do
    do rep = 1, nrep
      a2 = a2ref
      call cpu_time(t0)
      call man2(a2, mask, k0, m)
      call cpu_time(t1)
      tman2 = tman2 + (t1-t0)
    end do

    ! 1D integer
    tpacki = 0.0d0; tmani = 0.0d0
    do rep = 1, nrep
      ai = airef
      call cpu_time(t0)
      ai(1:m) = pack(ai(1:m), mask(1:m), ai(1:m))
      call cpu_time(t1)
      tpacki = tpacki + (t1-t0)
    end do
    do rep = 1, nrep
      ai = airef
      call cpu_time(t0)
      call mani(ai, mask, k0, m)
      call cpu_time(t1)
      tmani = tmani + (t1-t0)
    end do

    ! Verify manual == pack result
    a1 = a1ref
    call man1(a1, mask, 1, m)
    a1ref(1:m) = pack(a1ref(1:m), mask(1:m), a1ref(1:m))
    if (any(abs(a1(1:m-ndead) - a1ref(1:m-ndead)) > 0.0d0)) then
      print '(a)', 'ERROR: manual and pack results differ!'
      stop 1
    end if
    ai = airef
    call mani(ai, mask, k0, m)
    airef(1:m) = pack(airef(1:m), mask(1:m), airef(1:m))
    if (any(ai(1:m-ndead) /= airef(1:m-ndead))) then
      print '(a)', 'ERROR: manual(k0) and pack integer results differ!'
      stop 1
    end if

    print '(a, f10.3, a, f10.3, a, f6.2)', '  1D double:  pack ', tpack1*1d3/nrep, &
        ' ms   manual ', tman1*1d3/nrep, ' ms   speedup x', tpack1/tman1
    print '(a, f10.3, a, f10.3, a, f6.2)', '  2D double:  pack ', tpack2*1d3/nrep, &
        ' ms   manual ', tman2*1d3/nrep, ' ms   speedup x', tpack2/tman2
    print '(a, f10.3, a, f10.3, a, f6.2)', '  1D int:     pack ', tpacki*1d3/nrep, &
        ' ms   manual ', tmani*1d3/nrep, ' ms   speedup x', tpacki/tmani
    print '(a, i0)', '  (manual started at first dead index ', k0, ')'
    deallocate(a1, a1ref, a2, a2ref, ai, airef, mask)
  end do

contains

  subroutine make_mask(mask, m, ndead)
    logical, intent(out) :: mask(:)
    integer, intent(in)  :: m, ndead
    integer :: nd, idx
    double precision :: r
    mask = .true.
    nd = 0
    do
      if (nd >= ndead) exit
      call random_number(r)
      idx = 1 + int(r*m)
      if (mask(idx)) then
        mask(idx) = .false.
        nd = nd + 1
      end if
    end do
  end subroutine

  subroutine man1(A, mask, k, m)
    double precision, intent(inout) :: A(:)
    logical, intent(in) :: mask(:)
    integer, intent(in) :: k, m
    integer :: i, j
    j = k - 1
    do i = k, m
      if (mask(i)) then
        j = j + 1
        A(j) = A(i)
      end if
    end do
  end subroutine

  subroutine man2(A, mask, k, m)
    double precision, intent(inout) :: A(:,:)
    logical, intent(in) :: mask(:)
    integer, intent(in) :: k, m
    integer :: i, j
    j = k - 1
    do i = k, m
      if (mask(i)) then
        j = j + 1
        A(:, j) = A(:, i)
      end if
    end do
  end subroutine

  subroutine mani(A, mask, k, m)
    integer, intent(inout) :: A(:)
    logical, intent(in) :: mask(:)
    integer, intent(in) :: k, m
    integer :: i, j
    j = k - 1
    do i = k, m
      if (mask(i)) then
        j = j + 1
        A(j) = A(i)
      end if
    end do
  end subroutine

end program bench_compact
