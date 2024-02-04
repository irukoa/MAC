program Integer_Overflow_Happens

  use MAC, only: container_specifier
  implicit none

  type(container_specifier) :: specifier

  call specifier%specify([75, 75, 75, 75, 75])

end program Integer_Overflow_Happens
