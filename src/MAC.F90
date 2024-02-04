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

  pure subroutine set_al(self, value, at)
    class(container), intent(inout) :: self
    logical, intent(in) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 1) error stop &
      "MAC: Error #8: trying to set logical value to non logical container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%l_storage)) error stop &
      "MAC: Error #10: the logical container is not allocated."
    self%l_storage(memory_layout) = value
  end subroutine set_al

  pure subroutine set_ml(self, value, at)
    class(container), intent(inout) :: self
    logical, intent(in) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 1) error stop &
      "MAC: Error #8: trying to set logical value to non logical container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%l_storage)) error stop &
      "MAC: Error #10: the logical container is not allocated."
    self%l_storage(at) = value
  end subroutine set_ml

  pure subroutine set_gl(self, value)
    class(container), intent(inout) :: self
    logical, intent(in) :: value
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 1) error stop &
      "MAC: Error #8: trying to set logical value to non logical container."
    if (.not. allocated(self%l_storage)) error stop &
      "MAC: Error #10: the logical container is not allocated."
    self%l_storage = value
  end subroutine set_gl

  pure subroutine get_al(self, value, at)
    class(container), intent(inout) :: self
    logical, intent(out) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 1) error stop &
      "MAC: Error #8: trying to get logical value from non logical container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%l_storage)) error stop &
      "MAC: Error #10: the logical container is not allocated."
    value = self%l_storage(memory_layout)
  end subroutine get_al

  pure subroutine get_ml(self, value, at)
    class(container), intent(inout) :: self
    logical, intent(out) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 1) error stop &
      "MAC: Error #8: trying to get logical value from non logical container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%l_storage)) error stop &
      "MAC: Error #10: the logical container is not allocated."
    value = self%l_storage(at)
  end subroutine get_ml

  pure subroutine set_ai(self, value, at)
    class(container), intent(inout) :: self
    integer, intent(in) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 2) error stop &
      "MAC: Error #8: trying to set integer value to non integer container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%i_storage)) error stop &
      "MAC: Error #10: the integer container is not allocated."
    self%i_storage(memory_layout) = value
  end subroutine set_ai

  pure subroutine set_mi(self, value, at)
    class(container), intent(inout) :: self
    integer, intent(in) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 2) error stop &
      "MAC: Error #8: trying to set integer value to non integer container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%i_storage)) error stop &
      "MAC: Error #10: the integer container is not allocated."
    self%i_storage(at) = value
  end subroutine set_mi

  pure subroutine set_gi(self, value)
    class(container), intent(inout) :: self
    integer, intent(in) :: value
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 2) error stop &
      "MAC: Error #8: trying to set integer value to non integer container."
    if (.not. allocated(self%i_storage)) error stop &
      "MAC: Error #10: the integer container is not allocated."
    self%i_storage = value
  end subroutine set_gi

  pure subroutine get_ai(self, value, at)
    class(container), intent(inout) :: self
    integer, intent(out) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 2) error stop &
      "MAC: Error #8: trying to get integer value from non integer container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%i_storage)) error stop &
      "MAC: Error #10: the integer container is not allocated."
    value = self%i_storage(memory_layout)
  end subroutine get_ai

  pure subroutine get_mi(self, value, at)
    class(container), intent(inout) :: self
    integer, intent(out) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 2) error stop &
      "MAC: Error #8: trying to get integer value from non integer container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%i_storage)) error stop &
      "MAC: Error #10: the integer container is not allocated."
    value = self%i_storage(at)
  end subroutine get_mi

  pure subroutine set_ar(self, value, at)
    class(container), intent(inout) :: self
    real, intent(in) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 3) error stop &
      "MAC: Error #8: trying to set real value to non real container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%r_storage)) error stop &
      "MAC: Error #10: the real container is not allocated."
    self%r_storage(memory_layout) = value
  end subroutine set_ar

  pure subroutine set_mr(self, value, at)
    class(container), intent(inout) :: self
    real, intent(in) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 3) error stop &
      "MAC: Error #8: trying to set real value to non real container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%r_storage)) error stop &
      "MAC: Error #10: the real container is not allocated."
    self%r_storage(at) = value
  end subroutine set_mr

  pure subroutine set_gr(self, value)
    class(container), intent(inout) :: self
    real, intent(in) :: value
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 3) error stop &
      "MAC: Error #8: trying to set real value to non real container."
    if (.not. allocated(self%r_storage)) error stop &
      "MAC: Error #10: the real container is not allocated."
    self%r_storage = value
  end subroutine set_gr

  pure subroutine get_ar(self, value, at)
    class(container), intent(inout) :: self
    real, intent(out) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 3) error stop &
      "MAC: Error #8: trying to get real value from non real container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%r_storage)) error stop &
      "MAC: Error #10: the real container is not allocated."
    value = self%r_storage(memory_layout)
  end subroutine get_ar

  pure subroutine get_mr(self, value, at)
    class(container), intent(inout) :: self
    real, intent(out) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 3) error stop &
      "MAC: Error #8: trying to get real value from non real container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%r_storage)) error stop &
      "MAC: Error #10: the real container is not allocated."
    value = self%r_storage(at)
  end subroutine get_mr

  pure subroutine set_ac(self, value, at)
    class(container), intent(inout) :: self
    complex, intent(in) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 4) error stop &
      "MAC: Error #8: trying to set complex value to non complex container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%c_storage)) error stop &
      "MAC: Error #10: the complex container is not allocated."
    self%c_storage(memory_layout) = value
  end subroutine set_ac

  pure subroutine set_mc(self, value, at)
    class(container), intent(inout) :: self
    complex, intent(in) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 4) error stop &
      "MAC: Error #8: trying to set complex value to non complex container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%c_storage)) error stop &
      "MAC: Error #10: the complex container is not allocated."
    self%c_storage(at) = value
  end subroutine set_mc

  pure subroutine set_gc(self, value)
    class(container), intent(inout) :: self
    complex, intent(in) :: value
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 4) error stop &
      "MAC: Error #8: trying to set complex value to non complex container."
    if (.not. allocated(self%c_storage)) error stop &
      "MAC: Error #10: the complex container is not allocated."
    self%c_storage = value
  end subroutine set_gc

  pure subroutine get_ac(self, value, at)
    class(container), intent(inout) :: self
    complex, intent(out) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 4) error stop &
      "MAC: Error #8: trying to get complex value from non complex container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%c_storage)) error stop &
      "MAC: Error #10: the complex container is not allocated."
    value = self%c_storage(memory_layout)
  end subroutine get_ac

  pure subroutine get_mc(self, value, at)
    class(container), intent(inout) :: self
    complex, intent(out) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 4) error stop &
      "MAC: Error #8: trying to get complex value from non complex container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%c_storage)) error stop &
      "MAC: Error #10: the complex container is not allocated."
    value = self%c_storage(at)
  end subroutine get_mc

  pure subroutine set_ardp(self, value, at)
    class(container), intent(inout) :: self
    real(wp), intent(in) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 5) error stop &
      "MAC: Error #8: trying to set real(wp) value to non real(wp) container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%rdp_storage)) error stop &
      "MAC: Error #10: the real(wp) container is not allocated."
    self%rdp_storage(memory_layout) = value
  end subroutine set_ardp

  pure subroutine set_mrdp(self, value, at)
    class(container), intent(inout) :: self
    real(wp), intent(in) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 5) error stop &
      "MAC: Error #8: trying to set real(wp) value to non real(wp) container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%rdp_storage)) error stop &
      "MAC: Error #10: the real(wp) container is not allocated."
    self%rdp_storage(at) = value
  end subroutine set_mrdp

  pure subroutine set_grdp(self, value)
    class(container), intent(inout) :: self
    real(wp), intent(in) :: value
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 5) error stop &
      "MAC: Error #8: trying to set real(wp) value to non real(wp) container."
    if (.not. allocated(self%rdp_storage)) error stop &
      "MAC: Error #10: the real(wp) container is not allocated."
    self%rdp_storage = value
  end subroutine set_grdp

  pure subroutine get_ardp(self, value, at)
    class(container), intent(inout) :: self
    real(wp), intent(out) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 5) error stop &
      "MAC: Error #8: trying to get real(wp) value from non real(wp) container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%rdp_storage)) error stop &
      "MAC: Error #10: the real(wp) container is not allocated."
    value = self%rdp_storage(memory_layout)
  end subroutine get_ardp

  pure subroutine get_mrdp(self, value, at)
    class(container), intent(inout) :: self
    real(wp), intent(out) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 5) error stop &
      "MAC: Error #8: trying to get real(wp) value from non real(wp) container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%rdp_storage)) error stop &
      "MAC: Error #10: the real(wp) container is not allocated."
    value = self%rdp_storage(at)
  end subroutine get_mrdp

  pure subroutine set_acdp(self, value, at)
    class(container), intent(inout) :: self
    complex(wp), intent(in) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 6) error stop &
      "MAC: Error #8: trying to set complex(wp) value to non complex(wp) container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%cdp_storage)) error stop &
      "MAC: Error #10: the complex(wp) container is not allocated."
    self%cdp_storage(memory_layout) = value
  end subroutine set_acdp

  pure subroutine set_mcdp(self, value, at)
    class(container), intent(inout) :: self
    complex(wp), intent(in) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 6) error stop &
      "MAC: Error #8: trying to set complex(wp) value to non complex(wp) container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%cdp_storage)) error stop &
      "MAC: Error #10: the complex(wp) container is not allocated."
    self%cdp_storage(at) = value
  end subroutine set_mcdp

  pure subroutine set_gcdp(self, value)
    class(container), intent(inout) :: self
    complex(wp), intent(in) :: value
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 6) error stop &
      "MAC: Error #8: trying to set complex(wp) value to non complex(wp) container."
    if (.not. allocated(self%cdp_storage)) error stop &
      "MAC: Error #10: the complex(wp) container is not allocated."
    self%cdp_storage = value
  end subroutine set_gcdp

  pure subroutine get_acdp(self, value, at)
    class(container), intent(inout) :: self
    complex(wp), intent(out) :: value
    integer, intent(in) :: at(size(self%dimension_specifier))
    integer :: memory_layout
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 6) error stop &
      "MAC: Error #8: trying to get complex(wp) value from non complex(wp) container."
    memory_layout = self%ind(at)
    if ((memory_layout < 1) .or. (memory_layout > self%size())) error stop &
      "MAC: Error #9: specified array layout adress is out of bounds."
    if (.not. allocated(self%cdp_storage)) error stop &
      "MAC: Error #10: the complex(wp) container is not allocated."
    value = self%cdp_storage(memory_layout)
  end subroutine get_acdp

  pure subroutine get_mcdp(self, value, at)
    class(container), intent(inout) :: self
    complex(wp), intent(out) :: value
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 6) error stop &
      "MAC: Error #8: trying to get complex(wp) value from non complex(wp) container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%cdp_storage)) error stop &
      "MAC: Error #10: the complex(wp) container is not allocated."
    value = self%cdp_storage(at)
  end subroutine get_mcdp

end module MAC
