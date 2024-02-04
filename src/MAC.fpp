#:set types = ['logical', 'integer', 'real', 'complex', 'real(wp)', 'complex(wp)']
#:set numcode = [1, 2, 3, 4, 5, 6]
#:set storage_name = ['l_storage', 'i_storage', 'r_storage', 'c_storage', 'rdp_storage', 'cdp_storage']
#:set suffixes = ['l', 'i', 'r', 'c', 'rdp', 'cdp']
#:set lst = list(zip(types, numcode, storage_name, suffixes))
#:set interface_names = ['set', 'get']
module MAC

  use MAC_kinds, only: wp => dp

  implicit none
  private

  type, public :: container_specifier
    private
    integer, allocatable :: dimension_specifier(:)
    integer, allocatable :: lower_bounds(:)
    integer :: layout = 0
    logical :: specifier_initialized = .false.
  contains
    private
    procedure, public, pass(self) :: specify
    procedure, pass(self) :: memory_layout
    procedure, pass(self) :: array_layout
    generic, public :: ind => memory_layout, array_layout
    procedure, pass(self), public :: size => get_size
    procedure, pass(self) :: get_integer_property
    procedure, pass(self) :: get_integer_array_property
    generic, public :: get => get_integer_property, get_integer_array_property
    procedure, public, pass(self) :: partial_permutation
  end type

  type, extends(container_specifier), public :: container
    private
    integer :: container_type = 0
    logical :: container_initialized = .false.
    logical, allocatable, public :: l_storage(:)
    integer, allocatable, public :: i_storage(:)
    real, allocatable, public :: r_storage(:)
    complex, allocatable, public :: c_storage(:)
    real(wp), allocatable, public :: rdp_storage(:)
    complex(wp), allocatable, public :: cdp_storage(:)
  contains
    private
    procedure, public, pass(self) :: construct
    procedure, pass(self) :: set_al, set_ml, set_gl, get_al, get_ml
    procedure, pass(self) :: set_ai, set_mi, set_gi, get_ai, get_mi
    procedure, pass(self) :: set_ar, set_mr, set_gr, get_ar, get_mr
    procedure, pass(self) :: set_ac, set_mc, set_gc, get_ac, get_mc
    procedure, pass(self) :: set_ardp, set_mrdp, set_grdp, get_ardp, get_mrdp
    procedure, pass(self) :: set_acdp, set_mcdp, set_gcdp, get_acdp, get_mcdp
    generic, public :: set => set_al, set_ml, set_gl, &
      set_ai, set_mi, set_gi, &
      set_ar, set_mr, set_gr, &
      set_ac, set_mc, set_gc, &
      set_ardp, set_mrdp, set_grdp, &
      set_acdp, set_mcdp, set_gcdp
    generic, public :: get => get_al, get_ml, &
      get_ai, get_mi, &
      get_ar, get_mr, &
      get_ac, get_mc, &
      get_ardp, get_mrdp, &
      get_acdp, get_mcdp
  end type

contains

  pure subroutine partial_permutation(self, variables, dictionary)
    class(container_specifier), intent(in) :: self
    integer, intent(in) :: variables(:)
    integer, allocatable, intent(out) :: dictionary(:, :)

    character(len=1024) :: errormsg
    integer :: istat, i, n
    integer :: counter, reduction

    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    if (size(variables) > size(self%dimension_specifier)) error stop "MAC: Error #12: size of variables array is larger than rank."
    if (size(variables) == 0) error stop "MAC: Error #13: size of variables must be greater than 0."

    do i = 1, size(variables)
      if ((variables(i) < 1) .or. (variables(i) > size(self%dimension_specifier))) then
        write (errormsg, "(i20)") i
        errormsg = "MAC: Error #14: variables("//trim(adjustl(errormsg))//") does not reference a valid dimension index."
        error stop trim(errormsg)
      endif
    enddo

    n = 1
    do i = 1, size(variables)
      n = n*self%dimension_specifier(variables(i))
    enddo

    allocate (dictionary(n, size(self%dimension_specifier)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating permutation dictionary. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    dictionary = 1

    do i = 1, n

      associate (array_layout=>dictionary(i, :))
        select case (self%layout)
        case (0)
          reduction = i
          do counter = 1, size(variables) - 1
            array_layout(variables(counter)) = modulo(reduction, self%dimension_specifier(variables(counter)))
            if (array_layout(variables(counter)) == 0) &
              array_layout(variables(counter)) = self%dimension_specifier(variables(counter))
            reduction = int((reduction - array_layout(variables(counter)))/self%dimension_specifier(variables(counter))) + 1
          enddo
          array_layout(variables(counter)) = reduction
        case (1)
          reduction = i
          do counter = size(variables), 2, -1
            array_layout(variables(counter)) = modulo(reduction, self%dimension_specifier(variables(counter)))
            if (array_layout(variables(counter)) == 0) &
              array_layout(variables(counter)) = self%dimension_specifier(variables(counter))
            reduction = int((reduction - array_layout(variables(counter)))/self%dimension_specifier(variables(counter))) + 1
          enddo
          array_layout(variables(1)) = reduction
        end select
        array_layout = array_layout + self%lower_bounds - 1
      end associate

    enddo

  end subroutine

  subroutine specify(self, dimension_specifier, lower_bounds, layout)
    class(container_specifier), intent(out) :: self
    integer, intent(in) :: dimension_specifier(:)
    integer, intent(in), optional :: lower_bounds(:)
    character(len=*), intent(in), optional :: layout

    character(len=1024) :: errormsg
    integer :: istat, i

    if (lbound(dimension_specifier, dim=1) /= 1) error stop &
      "MAC: Error #11: the lower bound of dimension_specifier is not 1."

    if (present(lower_bounds)) then
      if (lbound(dimension_specifier, dim=1) /= 1) error stop &
        "MAC: Error #11: the lower bound of lower_bounds is not 1."
    endif

    if (present(lower_bounds)) then
      if (size(dimension_specifier) /= size(lower_bounds)) &
        error stop "MAC: Error #1: arrays dimension_specifier and lower_bounds are different sizes."
    endif

    do i = 1, size(dimension_specifier)
      if (dimension_specifier(i) < 1) then
        write (errormsg, "(i20)") i
        errormsg = "MAC: Error #2: dimension_specifier("//trim(adjustl(errormsg))//") is lower than 1."
        error stop trim(errormsg)
      endif
    enddo

    allocate (self%dimension_specifier(size(dimension_specifier)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating dimension_specifier. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    self%dimension_specifier = dimension_specifier

    safeguard_integer_overflow:block
      real(wp) :: tst
      integer :: i
      tst = 1.0_wp
      do i = 1, size(self%dimension_specifier)
        tst = tst*real(self%dimension_specifier(i), wp)
      enddo
      if (abs(tst - real(product(self%dimension_specifier), wp)) > epsilon(1.0_wp)) &
        error stop "MAC: Error #4: integer overflow. Requested array specifier too large."
    end block safeguard_integer_overflow

    allocate (self%lower_bounds(size(dimension_specifier)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating lower_bounds. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    if (present(lower_bounds)) then
      self%lower_bounds = lower_bounds
    else
      self%lower_bounds = 1
    endif

    if (present(layout)) then
      select case (layout)
      case ("column-major", "F", "col", "left")
        self%layout = 0
      case ("row-major", "C", "row", "right")
        self%layout = 1
      case default
        error stop "MAC: Error #5: specified layout not recognized."
      end select
    endif

    self%specifier_initialized = .true.

  end subroutine specify

  pure function array_layout(self, memory_layout)
    class(container_specifier), intent(in) :: self
    integer, intent(in) :: memory_layout
    integer :: array_layout(size(self%dimension_specifier))
    integer :: counter, reduction
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    array_layout = self%lower_bounds - 1
    if ((memory_layout < 1) .or. (memory_layout > product(self%dimension_specifier))) return
    array_layout = 0
    select case (self%layout)
    case (0)
      reduction = memory_layout
      do counter = 1, size(self%dimension_specifier) - 1
        array_layout(counter) = modulo(reduction, self%dimension_specifier(counter))
        if (array_layout(counter) == 0) array_layout(counter) = self%dimension_specifier(counter)
        reduction = int((reduction - array_layout(counter))/self%dimension_specifier(counter)) + 1
      enddo
      array_layout(size(self%dimension_specifier)) = reduction
    case (1)
      reduction = memory_layout
      do counter = size(self%dimension_specifier), 2, -1
        array_layout(counter) = modulo(reduction, self%dimension_specifier(counter))
        if (array_layout(counter) == 0) array_layout(counter) = self%dimension_specifier(counter)
        reduction = int((reduction - array_layout(counter))/self%dimension_specifier(counter)) + 1
      enddo
      array_layout(1) = reduction
    end select
    array_layout = array_layout + self%lower_bounds - 1
  end function array_layout

  pure integer function memory_layout(self, array_layout)
    class(container_specifier), intent(in) :: self
    integer, intent(in) :: array_layout(size(self%dimension_specifier))
    integer :: counter
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    memory_layout = 0
    do counter = 1, size(self%dimension_specifier)
      if ((array_layout(counter) < self%lower_bounds(counter)) .or. &
          (array_layout(counter) > self%dimension_specifier(counter) + self%lower_bounds(counter) - 1)) return
    enddo
    select case (self%layout)
    case (0)
      memory_layout = 1
      do counter = size(self%dimension_specifier), 1, -1
        memory_layout = (memory_layout - 1)*self%dimension_specifier(counter) + &
                        (array_layout(counter) - self%lower_bounds(counter) + 1)
      enddo
    case (1)
      memory_layout = 1
      do counter = 1, size(self%dimension_specifier)
        memory_layout = (memory_layout - 1)*self%dimension_specifier(counter) + &
                        (array_layout(counter) - self%lower_bounds(counter) + 1)
      enddo
    end select
  end function memory_layout

  pure elemental integer function get_size(self)
    class(container_specifier), intent(in) :: self
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    get_size = product(self%dimension_specifier)
  end function get_size

  pure elemental subroutine get_integer_property(self, property, val)
    class(container_specifier), intent(in) :: self
    character(len=*), intent(in) :: property
    integer, intent(out) :: val
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    select case (property)
    case ("size")
      val = product(self%dimension_specifier)
    case ("layout")
      val = self%layout
    case ("dimension", "rank")
      val = size(self%dimension_specifier)
    case default
      error stop "MAC: Error #7: requested property not recognized."
    end select
  end subroutine get_integer_property

  pure subroutine get_integer_array_property(self, property, val)
    class(container_specifier), intent(in) :: self
    character(len=*), intent(in) :: property
    integer, intent(out) :: val(size(self%dimension_specifier))
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    select case (property)
    case ("shape", "dimension_specifier")
      val = self%dimension_specifier
    case ("lower_bounds", "lb", "l_bounds")
      val = self%lower_bounds
    case ("upper_bounds", "ub", "u_bounds")
      val = self%lower_bounds + self%dimension_specifier - 1
    case default
      error stop "MAC: Error #7: requested property not recognized."
    end select
  end subroutine get_integer_array_property

  subroutine construct(self, container_type, dimension_specifier, lower_bounds, layout)
    class(container), intent(out) :: self
    character(len=*), intent(in) :: container_type
    integer, intent(in) :: dimension_specifier(:)
    integer, intent(in), optional :: lower_bounds(:)
    character(len=*), intent(in), optional :: layout

    character(len=1024) :: errormsg
    integer :: istat, i

    if (lbound(dimension_specifier, dim=1) /= 1) error stop &
      "MAC: Error #11: the lower bound of dimension_specifier is not 1."

    if (present(lower_bounds)) then
      if (lbound(dimension_specifier, dim=1) /= 1) error stop &
        "MAC: Error #11: the lower bound of lower_bounds is not 1."
    endif

    if (present(lower_bounds)) then
      if (size(dimension_specifier) /= size(lower_bounds)) &
        error stop "MAC: Error #1: array dimension_specifier and lower_bounds are different sizes."
    endif

    do i = 1, size(dimension_specifier)
      if (dimension_specifier(i) < 1) then
        write (errormsg, "(i20)") i
        errormsg = "MAC: Error #2: dimension_specifier("//trim(adjustl(errormsg))//") is lower than 1."
        error stop trim(errormsg)
      endif
    enddo

    allocate (self%dimension_specifier(size(dimension_specifier)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating dimension_specifier. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    self%dimension_specifier = dimension_specifier

    safeguard_integer_overflow:block
      real(wp) :: tst
      integer :: i
      tst = 1.0_wp
      do i = 1, size(self%dimension_specifier)
        tst = tst*real(self%dimension_specifier(i), wp)
      enddo
      if (abs(tst - real(product(self%dimension_specifier), wp)) > epsilon(1.0_wp)) &
        error stop "MAC: Error #4: integer overflow. Requested array specifier too large."
    end block safeguard_integer_overflow

    allocate (self%lower_bounds(size(dimension_specifier)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating lower_bounds. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    if (present(lower_bounds)) then
      self%lower_bounds = lower_bounds
    else
      self%lower_bounds = 1
    endif

    if (present(layout)) then
      select case (layout)
      case ("column-major", "F", "col", "left")
        self%layout = 0
      case ("row-major", "C", "row", "right")
        self%layout = 1
      case default
        error stop "MAC: Error #5: specified layout not recognized."
      end select
    endif

    self%specifier_initialized = .true.

    select case (container_type)
    case ("logical")
      self%container_type = 1
      allocate (self%l_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for logical container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%l_storage = .false.
    case ("integer")
      self%container_type = 2
      allocate (self%i_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for integer container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%i_storage = 0
    case ("real")
      self%container_type = 3
      allocate (self%r_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for real container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%r_storage = real(0.0)
    case ("complex")
      self%container_type = 4
      allocate (self%c_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for complex container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%c_storage = cmplx(0.0, 0.0)
    case ("real_dp")
      self%container_type = 5
      allocate (self%rdp_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for real container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%rdp_storage = real(0.0, wp)
    case ("complex_dp")
      self%container_type = 6
      allocate (self%cdp_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for complex container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%cdp_storage = cmplx(0.0, 0.0, wp)
    case default
      error stop "MAC: Error #5: requested container type not recognized."
    end select

    self%container_initialized = .true.

  end subroutine construct

#:for type, num, storage_name, suffix in lst
  pure subroutine set_a${suffix}$(self, value, at)
    class(container), intent(inout) :: self
    ${type}$, intent(in) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= ${num}$) error stop &
      "MAC: Error #8: trying to set ${type}$ value to non ${type}$ container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%${storage_name}$)) error stop &
      "MAC: Error #10: the ${type}$ container is not allocated."
    self%${storage_name}$(memory_layout) = value
  end subroutine set_a${suffix}$

  pure subroutine set_m${suffix}$(self, value, at)
    class(container), intent(inout) :: self
    ${type}$, intent(in) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= ${num}$) error stop &
      "MAC: Error #8: trying to set ${type}$ value to non ${type}$ container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%${storage_name}$)) error stop &
      "MAC: Error #10: the ${type}$ container is not allocated."
    self%${storage_name}$(at) = value
  end subroutine set_m${suffix}$

  pure subroutine set_g${suffix}$(self, value)
    class(container), intent(inout) :: self
    ${type}$, intent(in) :: value
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= ${num}$) error stop &
      "MAC: Error #8: trying to set ${type}$ value to non ${type}$ container."
    if (.not. allocated(self%${storage_name}$)) error stop &
      "MAC: Error #10: the ${type}$ container is not allocated."
    self%${storage_name}$ = value
  end subroutine set_g${suffix}$

  pure subroutine get_a${suffix}$(self, value, at)
    class(container), intent(inout) :: self
    ${type}$, intent(out) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= ${num}$) error stop &
      "MAC: Error #8: trying to get ${type}$ value from non ${type}$ container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%${storage_name}$)) error stop &
      "MAC: Error #10: the ${type}$ container is not allocated."
    value = self%${storage_name}$(memory_layout)
  end subroutine get_a${suffix}$

  pure subroutine get_m${suffix}$(self, value, at)
    class(container), intent(inout) :: self
    ${type}$, intent(out) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= ${num}$) error stop &
      "MAC: Error #8: trying to get ${type}$ value from non ${type}$ container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%${storage_name}$)) error stop &
      "MAC: Error #10: the ${type}$ container is not allocated."
    value = self%${storage_name}$(at)
  end subroutine get_m${suffix}$

#:endfor
end module MAC
