module Expected_Behaviour_Suite
  use MAC_kinds, only: wp => dp
  use MAC
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  public :: Collect_Expected_Behaviour_Tests

contains

  subroutine Collect_Expected_Behaviour_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Expected Behaviour", test_expected_bahaviour)]

  end subroutine Collect_Expected_Behaviour_Tests

  subroutine test_expected_bahaviour(error)
    type(error_type), allocatable, intent(out) :: error
    block
      type(container_specifier) :: a, b
      type(container) :: r
      integer :: i
      call a%specify(dimension_specifier=[5, 5, 5, 5, 5], &
                     lower_bounds=[-1, 3, 2, 7, 12], &
                     layout="F")
      call b%specify(dimension_specifier=[5, 5, 5, 5, 5], &
                     lower_bounds=[-1, 3, 2, 7, 12], &
                     layout="C")
      call r%construct(container_type="real_dp", &
                       dimension_specifier=[10, 10, 10, 10, 10], &
                       lower_bounds=[-1, 3, 2, 7, 12], &
                       layout="F")
      do i = 1, a%size()
        if (i /= a%ind(a%ind(i))) then
          allocate (error)
          return
        endif
        if (i /= b%ind(b%ind(i))) then
          allocate (error)
          return
        endif
      enddo
      do i = 1, r%size()
        call r%set(val=real(i, wp), at=r%ind(i))
      enddo
      if (abs((sum(r%rdp_storage)) - 100001.0_wp*50000.0_wp) > 50000.0_wp*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    end block
  end subroutine test_expected_bahaviour

end module Expected_Behaviour_Suite
