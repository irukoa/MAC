program example_basic

  !In this example we create a container_specifier (array handle),
  !and request some basic properties by using the API.

  use MAC, only: container_specifier
  use iso_fortran_env, only: output_unit

  implicit none

  type(container_specifier) :: a

  !Create a 2 x 3 x 4 3-dimensional array handle.
  call a%specify(dimension_specifier=[2, 3, 4])

  !Get its size.
  write (unit=output_unit, fmt=*) "Size = ", a%size(), "."

  !Get other properties and print them.
  write (unit=output_unit, fmt=*) "Rank = ", a%rank(), "."
  write (unit=output_unit, fmt=*) "Shape = ", a%shape(), "."
  write (unit=output_unit, fmt=*) "Lower Bounds = ", a%lbounds(), "."
  write (unit=output_unit, fmt=*) "Upper Bounds = ", a%ubounds(), "."
  write (unit=output_unit, fmt=*) "Upper Bounds = ", a%layout(), ", (0 = F, 1 = C)"

  !Now, we check what array layout adress corresponds to memory layout adresses 1, 5.
  !Since this is "F" layout, these should be [1, 1, 1] and [1, 3, 1] respectively.

  write (unit=output_unit, fmt=*) "Mem = 1, Arr = [", a%ind(1), "]."
  write (unit=output_unit, fmt=*) "Mem = 5, Arr = [", a%ind(5), "]."

  !Repeat for "C" layout. Now the array layout indices should be [1, 1, 1] and [1, 2, 1]
  !respectively.

  call a%specify(dimension_specifier=[2, 3, 4], layout="C")
  write (unit=output_unit, fmt=*) "Mem = 1, Arr = [", a%ind(1), "]."
  write (unit=output_unit, fmt=*) "Mem = 5, Arr = [", a%ind(5), "]."

end program example_basic
