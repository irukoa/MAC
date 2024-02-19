module Sections_Suite
  use, intrinsic :: iso_fortran_env, only: output_unit
  use MAC_kinds, only: wp => dp
  use MAC
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  public :: Collect_Sections_Tests

contains

  subroutine Collect_Sections_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Test array section (randomized):", test_section)]

  end subroutine Collect_Sections_Tests

  subroutine test_section(error)
    type(error_type), allocatable, intent(out) :: error

    type(container) :: a, b

    integer, parameter :: maxdim = 4, maxspec = 9, maxbounds = 200
    real(wp) :: random
    integer :: i, j, k, dim, ctype, tmp
    integer, allocatable :: spec(:), lbs(:), arr(:), sdims(:), sat(:), npds(:), seval(:)
    character(len=1) :: lyt
    logical :: cond

    logical :: lgl1, lgl2
    integer :: itg1, itg2
    real :: rls1, rls2
    complex :: cms1, cms2
    real(wp) :: rld1, rld2
    complex(wp) :: cmd1, cmd2

    call random_seed()
    call random_number(random)
    dim = nint(2.0_wp + real(maxdim - 1, wp)*random)

    !dim = 5 !DBG

    allocate (spec(dim), lbs(dim), arr(dim))
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

    !spec = [2, 3, 4, 5, 6] !DBG
    !lbs =  [-1, 0, 1, 1, -2]  !DBG
    !lyt = "F"!DBG
    !ctype = 3!DBG

    call a%construct(container_type=ctype, &
                     dimension_specifier=spec, &
                     lower_bounds=lbs, &
                     layout=lyt)

    select case (ctype)
    case (1)
      do i = 1, a%size()
        if (a%size() /= 1) then
          a%l_storage(i) = (real(i - 1, wp) > real(a%size(), wp)/2.0_wp)
        else
          a%l_storage(i) = .false.
        endif
      enddo
    case (2)
      do i = 1, a%size()
        a%i_storage(i) = i
      enddo
    case (3)
      do i = 1, a%size()
        a%r_storage(i) = real(i)
      enddo
    case (4)
      do i = 1, a%size()
        a%c_storage(i) = cmplx(real(i), 0.0)
      enddo
    case (5)
      do i = 1, a%size()
        a%rdp_storage(i) = real(i, wp)
      enddo
    case (6)
      do i = 1, a%size()
        a%cdp_storage(i) = cmplx(real(i, wp), 0.0_wp, wp)
      enddo
    end select

    call random_number(random)
    dim = nint(1.0_wp + real(size(spec) - 2, wp)*random)

    !dim = 5 !DBG

    allocate (sdims(dim), sat(size(spec) - dim), npds(size(spec) - dim), seval(dim))

    sdims = 0
    do i = 1, dim
2     continue
      call random_number(random)
      tmp = nint(1.0_wp + real(size(spec) - 1, wp)*random)
      if (any(sdims == tmp)) goto 2
      sdims(i) = tmp
      write (output_unit, fmt="(A, i2, A, i2, A)") &
      & "Dimension #", i, " chosen to reduce = ", sdims(i), "."
    enddo

    k = 1
    do i = 1, size(spec)
      cond = .false.
      if (any(sdims == i)) cond = .true.
      if (cond) then
        cycle
      else
        npds(k) = i
        k = k + 1
      endif
    enddo

    do i = 1, size(spec) - dim
      call random_number(random)
      sat(i) = nint(lbs(npds(i)) + real(spec(npds(i)) - 1, wp)*random)
      write (output_unit, fmt="(A, i2, A, i2, A, i3, A)") &
      & "Coordinate to evaluate reduction ", i, "/", size(spec) - dim, ". Value = ", sat(i), "."
    enddo

    b = a%section(dims=sdims, at=sat)

    do i = 1, dim
      call random_number(random)
      seval(i) = nint(lbs(sdims(i)) + real(spec(sdims(i)) - 1, wp)*random)
      write (output_unit, fmt="(A, i2, A, i2, A, i3, A)") &
      & "Coordinate to evaluate check ", i, "/", dim, ". Value = ", seval(i), "."
    enddo

    do i = 1, size(spec)
      if (any(sdims == i)) then
        j = findloc(sdims, i, 1)
        arr(i) = seval(j)
      else
        j = findloc(npds, i, 1)
        arr(i) = sat(j)
      endif
    enddo

    select case (b%cont_type())
    case (1)
      call b%get(lgl1, seval)
      call a%get(lgl2, arr)
      if (lgl1 .neqv. lgl2) then
        allocate (error)
        return
      endif
    case (2)
      call b%get(itg1, seval)
      call a%get(itg2, arr)
      if (itg1 /= itg2) then
        allocate (error)
        return
      endif
    case (3)
      call b%get(rls1, seval)
      call a%get(rls2, arr)
      if (abs(rls1 - rls2) > 1.0E1*epsilon(1.0)) then
        allocate (error)
        return
      endif
    case (4)
      call b%get(cms1, seval)
      call a%get(cms2, arr)
      if (abs(abs(cms1) - abs(cms2)) > 1.0E1*epsilon(1.0)) then
        allocate (error)
        return
      endif
    case (5)
      call b%get(rld1, seval)
      call a%get(rld2, arr)
      if (abs(rld1 - rld2) > 1.0E1_wp*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    case (6)
      call b%get(cmd1, seval)
      call a%get(cmd2, arr)
      if (abs(abs(cmd1) - abs(cmd2)) > 1.0E1_wp*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    end select

  end subroutine test_section

end module Sections_Suite
