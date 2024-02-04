module Randomized_Suite
  use, intrinsic :: iso_fortran_env, only: output_unit
  use MAC_kinds, only: wp => dp
  use MAC
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  public :: Collect_Randomized_Tests

contains

  subroutine Collect_Randomized_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Randomized:", test_randomized)]

  end subroutine Collect_Randomized_Tests

  subroutine test_randomized(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: maxdim = 10, maxspec = 9, maxbounds = 200
    block
      type(container) :: r
      real(wp) :: random, check
      integer :: i, dim
      integer, allocatable :: spec(:), lbs(:)
      character(len=1) :: lyt

      call random_init(.true., .true.)
      call random_seed()
      call random_number(random)
      dim = nint(1.0_wp + real(maxdim - 1, wp)*random)
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
      call r%construct(container_type="real_dp", &
                       dimension_specifier=spec, &
                       lower_bounds=lbs, &
                       layout=lyt)
      do i = 1, r%size()
        call r%set(value=log(real(i, wp)), at=r%ind(r%ind(i)))
      enddo
      call random_number(random)
      i = nint(1.0_wp + r%size()*random)
      call r%get(value=check, at=r%ind(i))
      if (abs(real(i, wp) - exp(check)) > 1.0E10_wp*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
      deallocate (spec, lbs)
    end block
  end subroutine test_randomized

end module Randomized_Suite
