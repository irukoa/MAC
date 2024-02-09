program example_partial_permutation

  !In this example, we show the partial permutation utility.
  !This feature allows to permute only on a subset of dimensions
  !and get the array layout adress of said permutations.

  use MAC, only: container_specifier
  use iso_fortran_env, only: output_unit

  implicit none

  type(container_specifier) :: a
  integer :: i
  integer, allocatable :: dict(:, :)

  !Create a 2 x 3 x 4 3-dimensional array handle.
  call a%specify(dimension_specifier=[2, 3, 4], layout="F")

  !Consider the permutation of only dimensions #1 and #2.
  !The permutations are: (1, 1, 1), (2, 1, 1), (1, 2, 1)...
  dict = a%partial_permutation(dims=[1, 2])
  do i = 1, size(dict(:, 1))
    write (unit=output_unit, fmt=*) "Perm = ", i, ". Arr = ", dict(i, :), ", Mem = ", a%ind(dict(i, :)), "."
  enddo
  !so dimension #3 is not permuted over.

  !Notice that the order of the specified dimensions to permute over matters,
  !`a` has "F" layout, so the leftmost dimension label specified in `dims` is what gets
  !permuted first. For example,
  dict = a%partial_permutation(dims=[2, 1])
  do i = 1, size(dict(:, 1))
    write (unit=output_unit, fmt=*) "Perm = ", i, ". Arr = ", dict(i, :), ", Mem = ", a%ind(dict(i, :)), "."
  enddo
  !permutes first dimension label #2 and then dimension label #1.

  deallocate (dict)

end program example_partial_permutation
