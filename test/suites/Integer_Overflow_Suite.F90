module Integer_Overflow_Suite
  use MAC
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  public :: Collect_Integer_Overflow_Tests

contains

  subroutine Collect_Integer_Overflow_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Integer Overflow test", test_integer_overflow_does_not_happen), &
                 new_unittest("Integer Overflow test", test_integer_overflow_happens, should_fail=.true.)]

  end subroutine Collect_Integer_Overflow_Tests

  subroutine test_integer_overflow_does_not_happen(error)
    type(error_type), allocatable, intent(out) :: error
    type(container_specifier) :: a
    integer :: i
    call execute_command_line("fpm run Integer_Overflow_Does_Not_Happen > out.log 2>&1", wait=.true., exitstat=i)
    if (i == 1) allocate (error)
    if (allocated(error)) return
  end subroutine test_integer_overflow_does_not_happen

  subroutine test_integer_overflow_happens(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    call execute_command_line("fpm run Integer_Overflow_Happens > out.log 2>&1", wait=.true., exitstat=i)
    if (i == 1) allocate (error)
    if (allocated(error)) return
  end subroutine test_integer_overflow_happens

end module Integer_Overflow_Suite
