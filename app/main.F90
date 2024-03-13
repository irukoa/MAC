program main

  use MAC, only: container_specifier
  use iso_fortran_env, only: output_unit

  implicit none

  character(len=10) :: ver = "1.0.0"

  write (unit=output_unit, fmt="(A)") ""
  write (unit=output_unit, fmt="(A)") " .----------------.  .----------------.  .----------------.   "
  write (unit=output_unit, fmt="(A)") " | .--------------. || .--------------. || .--------------. | "
  write (unit=output_unit, fmt="(A)") " | | ____    ____ | || |      __      | || |     ______   | | "
  write (unit=output_unit, fmt="(A)") " | ||_   \  /   _|| || |     /  \     | || |   .' ___  |  | | "
  write (unit=output_unit, fmt="(A)") " | |  |   \/   |  | || |    / /\ \    | || |  / .'   \_|  | | "
  write (unit=output_unit, fmt="(A)") " | |  | |\  /| |  | || |   / ____ \   | || |  | |         | | "
  write (unit=output_unit, fmt="(A)") " | | _| |_\/_| |_ | || | _/ /    \ \_ | || |  \ `.___.'\  | | "
  write (unit=output_unit, fmt="(A)") " | ||_____||_____|| || ||____|  |____|| || |   `._____.'  | | "
  write (unit=output_unit, fmt="(A)") " | |              | || |              | || |              | | "
  write (unit=output_unit, fmt="(A)") " | '--------------' || '--------------' || '--------------' | "
  write (unit=output_unit, fmt="(A)") "  '----------------'  '----------------'  '----------------'  "
  write (unit=output_unit, fmt="(A)") ""
  write (unit=output_unit, fmt="(A)") " Multidimensional Array Containers (MAC) v"//trim(adjustl(ver))//" built."

end program main
