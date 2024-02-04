program Integer_Overflow_Does_Not_Happen

  use MAC, only: container_specifier
  implicit none

  type(container_specifier) :: specifier

  call specifier%specify([75, 75, 75, 75, 67])

end program Integer_Overflow_Does_Not_Happen
