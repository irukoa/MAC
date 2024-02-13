program example_container

  !In this example we create a container,
  !and reshape it to regular arrays.

  use MAC, only: container
  use iso_fortran_env, only: output_unit

  implicit none

  type(container) :: a
  real :: b(2, 3, 4)
  integer :: i1, i2, i3

  !Create a 2 x 3 x 4 3-dimensional real array.
  call a%construct(container_type="real", dimension_specifier=[2, 3, 4])

  do i3 = 1, 4
    do i2 = 1, 3
      do i1 = 1, 2
        call a%set(val=real(i1*i2*i3), at=[i1, i2, i3])
      enddo
    enddo
  enddo

  b = reshape(a%r_storage, [2, 3, 4])

  do i3 = 1, 4
    do i2 = 1, 3
      do i1 = 1, 2
        write (unit=output_unit, fmt=*) b(i1, i2, i3)
      enddo
    enddo
  enddo

end program example_container
