module Reduction_Suite
  use, intrinsic :: iso_fortran_env, only: output_unit
  use MAC_kinds, only: wp => dp
  use MAC
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  integer, parameter :: tol = 1.0E2_wp

  public :: Collect_Reduction_Tests

contains

  subroutine Collect_Reduction_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Reduce a single dimension (randomized):", test_dim_reduction), &
                 new_unittest("Reduce all dimensions (randomized):", test_all_reduction)]

  end subroutine Collect_Reduction_Tests

  subroutine test_dim_reduction(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: maxdim = 4, maxspec = 9, maxbounds = 200
    real(wp) :: random, check
    integer :: i, dim, ctype, cdim
    integer, allocatable :: spec(:), lbs(:)
    character(len=1) :: lyt

    type(container) :: a
    real(wp) :: sm

    logical :: lgl
    integer :: itg
    real :: rls
    complex :: cms
    real(wp) :: rld
    complex(wp) :: cmd

    call random_seed()
    call random_number(random)
    dim = nint(2.0_wp + real(maxdim - 1, wp)*random)

    !dim = 5 !DBG

    allocate (spec(dim), lbs(dim))
    do i = 1, dim
      call random_number(random)
      spec(i) = nint(1.0_wp + real(maxspec - 1, wp)*random)
      call random_number(random)
      lbs(i) = nint(-100.0_wp + real(maxbounds - 1, wp)*random)
      write (output_unit, fmt="(A, i2, A, i2, A, i1, A, i3, A)") &
      & "Dimension, ", i, "/", dim, ": size = ", spec(i), ", lbound = ", lbs(i), "."
    enddo
    call random_number(random)
    if ((random - 0.5_wp) <= 0.0_wp) then
      lyt = "F"
    else
      lyt = "C"
    endif
    write (output_unit, fmt="(A)") "Layout = "//lyt//"."
    call random_number(random)
    ctype = nint(1.0_wp + 5.0_wp*random)
    write (output_unit, fmt="(A, i0, A)") "Container Type = ", ctype, "."

    !spec = 9 !DBG
    !lbs = 1  !DBG
    !lyt = "F"!DBG
    !ctype = 1!DBG

    call a%construct(container_type=ctype, &
                     dimension_specifier=spec, &
                     lower_bounds=lbs, &
                     layout=lyt)

    call random_number(random)
    cdim = nint(1.0_wp + real(dim - 1, wp)*random)
    write (output_unit, fmt="(A, i0, A)") "Dimension chosen to reduce = ", cdim, "."

    sm = real(spec(cdim), wp)

    select case (ctype)
    case (1)
      do i = 1, a%size()
        if (a%size() /= 1) then
          a%l_storage(i) = (real(i - 1, wp) > real(a%size(), wp)/2.0_wp)
        else
          a%l_storage(i) = .false.
        endif
      enddo
      a = a%reduce(sl, cdim)
      lgl = a%l_storage(1)
      if (lgl) then
        allocate (error)
        return
      endif
    case (2)
      do i = 1, a%size()
        a%i_storage(i) = 1
      enddo
      a = a%reduce(si, cdim)
      itg = a%i_storage(1)
      if (abs(real(itg, wp) - sm) > tol*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    case (3)
      do i = 1, a%size()
        a%r_storage(i) = 1.0
      enddo
      a = a%reduce(sr, cdim)
      rls = a%r_storage(1)
      if (abs(real(rls) - real(sm)) > tol*epsilon(1.0)) then
        allocate (error)
        return
      endif
    case (4)
      do i = 1, a%size()
        a%c_storage(i) = cmplx(1.0, 0.0)
      enddo
      a = a%reduce(sc, cdim)
      cms = a%c_storage(1)
      if (abs(real(cms) - real(sm)) > tol*epsilon(1.0)) then
        allocate (error)
        return
      endif
    case (5)
      do i = 1, a%size()
        a%rdp_storage(i) = 1.0_wp
      enddo
      a = a%reduce(srdp, cdim)
      rld = a%rdp_storage(1)
      if (abs(rld - sm) > tol*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    case (6)
      do i = 1, a%size()
        a%cdp_storage(i) = cmplx(1.0_wp, 0.0_wp, wp)
      enddo
      a = a%reduce(scdp, cdim)
      cmd = a%cdp_storage(1)
      if (abs(real(cmd, wp) - sm) > tol*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    end select

    deallocate (spec, lbs)

  contains
    function sl(array)
      logical, intent(in) :: array(:)
      logical :: sl
      sl = array(1)
    end function sl
    function si(array)
      integer, intent(in) :: array(:)
      integer :: si
      si = sum(array)
    end function si
    function sr(array)
      real, intent(in) :: array(:)
      real :: sr
      sr = sum(array)
    end function sr
    function sc(array)
      complex, intent(in) :: array(:)
      complex :: sc
      sc = sum(array)
    end function sc
    function srdp(array)
      real(wp), intent(in) :: array(:)
      real(wp) :: srdp
      srdp = sum(array)
    end function srdp
    function scdp(array)
      complex(wp), intent(in) :: array(:)
      complex(wp) :: scdp
      scdp = sum(array)
    end function scdp
  end subroutine test_dim_reduction

  subroutine test_all_reduction(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: maxdim = 5, maxspec = 9, maxbounds = 200
    real(wp) :: random, check
    integer :: i, dim, ctype
    integer, allocatable :: spec(:), lbs(:)
    character(len=1) :: lyt

    type(container) :: a
    real(wp) :: sm

    logical :: lgl
    integer :: itg
    real :: rls
    complex :: cms
    real(wp) :: rld
    complex(wp) :: cmd

    call random_seed()
    call random_number(random)
    dim = nint(1.0_wp + real(maxdim - 1, wp)*random)

    !dim = 5 !DBG

    allocate (spec(dim), lbs(dim))
    do i = 1, dim
      call random_number(random)
      spec(i) = nint(1.0_wp + real(maxspec - 1, wp)*random)
      call random_number(random)
      lbs(i) = nint(-100.0_wp + real(maxbounds - 1, wp)*random)
      write (output_unit, fmt="(A, i2, A, i2, A, i1, A, i3, A)") &
      & "Dimension, ", i, "/", dim, ": size = ", spec(i), ", lbound = ", lbs(i), "."
    enddo
    call random_number(random)
    if ((random - 0.5_wp) <= 0.0_wp) then
      lyt = "F"
    else
      lyt = "C"
    endif
    write (output_unit, fmt="(A)") "Layout = "//lyt//"."
    call random_number(random)
    ctype = nint(1.0_wp + 5.0_wp*random)
    write (output_unit, fmt="(A, i0, A)") "Container Type = ", ctype, "."

    !spec = 9 !DBG
    !lbs = 1  !DBG
    !lyt = "F"!DBG
    !ctype = 1!DBG

    call a%construct(container_type=ctype, &
                     dimension_specifier=spec, &
                     lower_bounds=lbs, &
                     layout=lyt)

    sm = 0.0_wp
    do i = 1, a%size()
      sm = sm + real(i, wp)
    enddo

    select case (ctype)
    case (1)
      do i = 1, a%size()
        if (a%size() /= 1) then
          a%l_storage(i) = (real(i - 1, wp) > real(a%size(), wp)/2.0_wp)
        else
          a%l_storage(i) = .false.
        endif
      enddo
      lgl = a%reduce(sl)
      if (lgl) then
        allocate (error)
        return
      endif
    case (2)
      do i = 1, a%size()
        a%i_storage(i) = i
      enddo
      itg = a%reduce(si)
      if (abs(real(itg, wp) - sm) > tol*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    case (3)
      do i = 1, a%size()
        a%r_storage(i) = real(i)
      enddo
      rls = a%reduce(sr)
      if (abs(real(rls) - real(sm)) > tol*epsilon(1.0)) then
        allocate (error)
        return
      endif
    case (4)
      do i = 1, a%size()
        a%c_storage(i) = cmplx(real(i), 0.0)
      enddo
      cms = a%reduce(sc)
      if (abs(real(cms) - real(sm)) > tol*epsilon(1.0)) then
        allocate (error)
        return
      endif
    case (5)
      do i = 1, a%size()
        a%rdp_storage(i) = real(i, wp)
      enddo
      rld = a%reduce(srdp)
      if (abs(rld - sm) > tol*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    case (6)
      do i = 1, a%size()
        a%cdp_storage(i) = cmplx(real(i, wp), 0.0_wp, wp)
      enddo
      cmd = a%reduce(scdp)
      if (abs(real(cmd, wp) - sm) > tol*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    end select

    deallocate (spec, lbs)

  contains
    function sl(array)
      logical, intent(in) :: array(:)
      logical :: sl
      sl = array(1)
    end function sl
    function si(array)
      integer, intent(in) :: array(:)
      integer :: si
      si = sum(array)
    end function si
    function sr(array)
      real, intent(in) :: array(:)
      real :: sr
      sr = sum(array)
    end function sr
    function sc(array)
      complex, intent(in) :: array(:)
      complex :: sc
      sc = sum(array)
    end function sc
    function srdp(array)
      real(wp), intent(in) :: array(:)
      real(wp) :: srdp
      srdp = sum(array)
    end function srdp
    function scdp(array)
      complex(wp), intent(in) :: array(:)
      complex(wp) :: scdp
      scdp = sum(array)
    end function scdp
  end subroutine test_all_reduction

end module Reduction_Suite
