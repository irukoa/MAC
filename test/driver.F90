program driver
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type, &
  & select_suite, run_selected, get_argument
  !Tests:
  use Integer_Overflow_Suite, only: Collect_Integer_Overflow_Tests
  use Expected_Behaviour_Suite, only: Collect_Expected_Behaviour_Tests
  use Randomized_Suite, only: Collect_Randomized_Tests
  use Reduction_Suite, only: Collect_Reduction_Tests

  implicit none
  integer :: stat, is
  character(len=:), allocatable :: suite_name, test_name
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [new_testsuite("Integer_Overflow", Collect_Integer_Overflow_Tests), &
                new_testsuite("Expected Behaviour", Collect_Expected_Behaviour_Tests), &
                new_testsuite("Randomized", Collect_Randomized_Tests), &
                new_testsuite("Reduction", Collect_Reduction_Tests)]

  call get_argument(1, suite_name)
  call get_argument(2, test_name)

  if (allocated(suite_name)) then
    is = select_suite(testsuites, suite_name)
    if (is > 0 .and. is <= size(testsuites)) then
      if (allocated(test_name)) then
        write (error_unit, fmt) "Suite:", testsuites(is)%name
        call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
        if (stat < 0) then
          error stop 1
        end if
      else
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end if
    else
      write (error_unit, fmt) "Available testsuites"
      do is = 1, size(testsuites)
        write (error_unit, fmt) "-", testsuites(is)%name
      end do
      error stop 1
    end if
  else
    do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do
  end if

  if (stat > 0) then
    write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if

end program driver
