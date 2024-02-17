module MAC

  use MAC_kinds, only: wp => dp

  implicit none
  private

  type, public :: container_specifier
    private
    integer, allocatable :: dimension_specifier(:)
    integer, allocatable :: lower_bounds(:)
    integer :: lyt = 0
    logical :: specifier_initialized = .false.
  contains
    private
    procedure, pass(self) :: spec_copy
    generic, public :: assignment(=) => spec_copy
    procedure, pass(self) :: specify_char
    procedure, pass(self) :: specify_int
    generic, public :: specify => specify_char, specify_int
    procedure, pass(self) :: array_layout
    procedure, pass(self) :: memory_layout
    generic, public :: ind => memory_layout, array_layout
    procedure, pass(self), public :: size => get_size
    procedure, pass(self), public :: rank => get_rank
    procedure, pass(self), public :: shape => get_shape
    procedure, pass(self), public :: lbounds => get_lbounds
    procedure, pass(self), public :: ubounds => get_ubounds
    procedure, pass(self), public :: layout => get_layout
    procedure, pass(self), public :: spec_init => get_is_specifier_initailized
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
    procedure, pass(self) :: cont_copy
    generic, public :: assignment(=) => cont_copy
    procedure, pass(self) :: construct_ctc_lc, construct_cti_lc, &
      construct_ctc_li, construct_cti_li
    generic, public :: construct => construct_ctc_lc, construct_cti_lc, &
      construct_ctc_li, construct_cti_li
    procedure, pass(self), public :: cont_type => get_container_type
    procedure, pass(self), public :: cont_init => get_is_container_initailized
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
    procedure, pass(self) :: rd_l, rd_i, rd_r, rd_c, rd_rdp, rd_cdp, &
      ra_l, ra_i, ra_r, ra_c, ra_rdp, ra_cdp
    generic, public :: reduce => rd_l, rd_i, rd_r, rd_c, rd_rdp, rd_cdp, &
      ra_l, ra_i, ra_r, ra_c, ra_rdp, ra_cdp
  end type

contains

  subroutine spec_copy(self, from)
    class(container_specifier), intent(out) :: self
    type(container_specifier), intent(in) :: from
    if (.not. from%spec_init()) error stop "MAC: Error #6: container specifier not initalized."
    call self%specify(from%shape(), from%lbounds(), from%layout())
  end subroutine spec_copy

  subroutine specify_char(self, dimension_specifier, lower_bounds, layout)
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
        self%lyt = 0
      case ("row-major", "C", "row", "right")
        self%lyt = 1
      case default
        error stop "MAC: Error #5: specified layout not recognized."
      end select
    endif

    self%specifier_initialized = .true.

  end subroutine specify_char

  subroutine specify_int(self, dimension_specifier, lower_bounds, layout)
    class(container_specifier), intent(out) :: self
    integer, intent(in) :: dimension_specifier(:)
    integer, intent(in), optional :: lower_bounds(:)
    integer, intent(in) :: layout

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

    select case (layout)
    case (0)
      self%lyt = 0
    case (1)
      self%lyt = 1
    case default
      error stop "MAC: Error #5: specified layout not recognized."
    end select

    self%specifier_initialized = .true.

  end subroutine specify_int

  pure function array_layout(self, memory_layout)
    class(container_specifier), intent(in) :: self
    integer, intent(in) :: memory_layout
    integer :: array_layout(size(self%dimension_specifier))
    integer :: counter, reduction
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    array_layout = self%lower_bounds - 1
    if ((memory_layout < 1) .or. (memory_layout > product(self%dimension_specifier))) return
    array_layout = 0
    select case (self%lyt)
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
    select case (self%lyt)
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

  pure elemental integer function get_rank(self)
    class(container_specifier), intent(in) :: self
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    get_rank = size(self%dimension_specifier)
  end function get_rank

  pure function get_shape(self)
    class(container_specifier), intent(in) :: self
    integer :: get_shape(size(self%dimension_specifier))
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    get_shape = self%dimension_specifier
  end function get_shape

  pure function get_lbounds(self)
    class(container_specifier), intent(in) :: self
    integer :: get_lbounds(size(self%dimension_specifier))
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    get_lbounds = self%lower_bounds
  end function get_lbounds

  pure function get_ubounds(self)
    class(container_specifier), intent(in) :: self
    integer :: get_ubounds(size(self%dimension_specifier))
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    get_ubounds = self%lower_bounds + self%dimension_specifier - 1
  end function get_ubounds

  pure elemental integer function get_layout(self)
    class(container_specifier), intent(in) :: self
    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    get_layout = self%lyt
  end function get_layout

  pure elemental logical function get_is_specifier_initailized(self)
    class(container_specifier), intent(in) :: self
    get_is_specifier_initailized = self%specifier_initialized
  end function get_is_specifier_initailized

  pure function partial_permutation(self, dims) result(dictionary)
    class(container_specifier), intent(in) :: self
    integer, intent(in) :: dims(:)
    integer, allocatable :: dictionary(:, :)

    character(len=1024) :: errormsg
    integer :: istat, i, n
    integer :: counter, reduction

    if (.not. self%specifier_initialized) error stop "MAC: Error #6: container specifier not initalized."
    if (size(dims) > size(self%dimension_specifier)) error stop "MAC: Error #12: size of dims array is larger than rank."
    if (size(dims) == 0) error stop "MAC: Error #7: size of dims must be greater than 0."

    do i = 1, size(dims)
      if ((dims(i) < 1) .or. (dims(i) > size(self%dimension_specifier))) then
        write (errormsg, "(i20)") i
        errormsg = "MAC: Error #14: dims("//trim(adjustl(errormsg))//") does not reference a valid dimension label."
        error stop trim(errormsg)
      endif
    enddo

    n = 1
    do i = 1, size(dims)
      n = n*self%dimension_specifier(dims(i))
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
        select case (self%lyt)
        case (0)
          reduction = i
          do counter = 1, size(dims) - 1
            array_layout(dims(counter)) = modulo(reduction, self%dimension_specifier(dims(counter)))
            if (array_layout(dims(counter)) == 0) &
              array_layout(dims(counter)) = self%dimension_specifier(dims(counter))
            reduction = int((reduction - array_layout(dims(counter)))/self%dimension_specifier(dims(counter))) + 1
          enddo
          array_layout(dims(counter)) = reduction
        case (1)
          reduction = i
          do counter = size(dims), 2, -1
            array_layout(dims(counter)) = modulo(reduction, self%dimension_specifier(dims(counter)))
            if (array_layout(dims(counter)) == 0) &
              array_layout(dims(counter)) = self%dimension_specifier(dims(counter))
            reduction = int((reduction - array_layout(dims(counter)))/self%dimension_specifier(dims(counter))) + 1
          enddo
          array_layout(dims(1)) = reduction
        end select
        array_layout = array_layout + self%lower_bounds - 1
      end associate

    enddo

  end function partial_permutation

  subroutine cont_copy(self, from)
    class(container), intent(out) :: self
    class(container), intent(in) :: from
    if (.not. from%cont_init()) error stop "MAC: Error #6: container not initalized."
    call self%construct(from%cont_type(), from%shape(), from%lbounds(), from%layout())
    select case (from%cont_type())
    case (1)
      self%l_storage = from%l_storage
    case (2)
      self%i_storage = from%i_storage
    case (3)
      self%r_storage = from%r_storage
    case (4)
      self%c_storage = from%c_storage
    case (5)
      self%rdp_storage = from%rdp_storage
    case (6)
      self%cdp_storage = from%cdp_storage
    end select
  end subroutine cont_copy

  subroutine construct_ctc_lc(self, container_type, dimension_specifier, lower_bounds, layout)
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
        self%lyt = 0
      case ("row-major", "C", "row", "right")
        self%lyt = 1
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

  end subroutine construct_ctc_lc

  subroutine construct_cti_lc(self, container_type, dimension_specifier, lower_bounds, layout)
    class(container), intent(out) :: self
    integer, intent(in) :: container_type
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
        self%lyt = 0
      case ("row-major", "C", "row", "right")
        self%lyt = 1
      case default
        error stop "MAC: Error #5: specified layout not recognized."
      end select
    endif

    self%specifier_initialized = .true.

    select case (container_type)
    case (1)
      self%container_type = 1
      allocate (self%l_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for logical container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%l_storage = .false.
    case (2)
      self%container_type = 2
      allocate (self%i_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for integer container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%i_storage = 0
    case (3)
      self%container_type = 3
      allocate (self%r_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for real container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%r_storage = real(0.0)
    case (4)
      self%container_type = 4
      allocate (self%c_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for complex container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%c_storage = cmplx(0.0, 0.0)
    case (5)
      self%container_type = 5
      allocate (self%rdp_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for real container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%rdp_storage = real(0.0, wp)
    case (6)
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

  end subroutine construct_cti_lc

  subroutine construct_ctc_li(self, container_type, dimension_specifier, lower_bounds, layout)
    class(container), intent(out) :: self
    character(len=*), intent(in) :: container_type
    integer, intent(in) :: dimension_specifier(:)
    integer, intent(in), optional :: lower_bounds(:)
    integer, intent(in) :: layout

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

    select case (layout)
    case (0)
      self%lyt = 0
    case (1)
      self%lyt = 1
    case default
      error stop "MAC: Error #5: specified layout not recognized."
    end select

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

  end subroutine construct_ctc_li

  subroutine construct_cti_li(self, container_type, dimension_specifier, lower_bounds, layout)
    class(container), intent(out) :: self
    integer, intent(in) :: container_type
    integer, intent(in) :: dimension_specifier(:)
    integer, intent(in), optional :: lower_bounds(:)
    integer, intent(in) :: layout

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

    select case (layout)
    case (0)
      self%lyt = 0
    case (1)
      self%lyt = 1
    case default
      error stop "MAC: Error #5: specified layout not recognized."
    end select

    self%specifier_initialized = .true.

    select case (container_type)
    case (1)
      self%container_type = 1
      allocate (self%l_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for logical container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%l_storage = .false.
    case (2)
      self%container_type = 2
      allocate (self%i_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for integer container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%i_storage = 0
    case (3)
      self%container_type = 3
      allocate (self%r_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for real container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%r_storage = real(0.0)
    case (4)
      self%container_type = 4
      allocate (self%c_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for complex container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%c_storage = cmplx(0.0, 0.0)
    case (5)
      self%container_type = 5
      allocate (self%rdp_storage(product(self%dimension_specifier)), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "MAC: Error #3: failure allocating storage for real container. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      self%rdp_storage = real(0.0, wp)
    case (6)
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

  end subroutine construct_cti_li

  pure elemental integer function get_container_type(self)
    class(container), intent(in) :: self
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    get_container_type = self%container_type
  end function get_container_type

  pure elemental logical function get_is_container_initailized(self)
    class(container), intent(in) :: self
    get_is_container_initailized = self%container_initialized
  end function get_is_container_initailized

  pure subroutine set_al(self, val, at)
    class(container), intent(inout) :: self
    logical, intent(in) :: val
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
    self%l_storage(memory_layout) = val
  end subroutine set_al

  pure subroutine set_ml(self, val, at)
    class(container), intent(inout) :: self
    logical, intent(in) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 1) error stop &
      "MAC: Error #8: trying to set logical value to non logical container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%l_storage)) error stop &
      "MAC: Error #10: the logical container is not allocated."
    self%l_storage(at) = val
  end subroutine set_ml

  pure subroutine set_gl(self, val)
    class(container), intent(inout) :: self
    logical, intent(in) :: val
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 1) error stop &
      "MAC: Error #8: trying to set logical value to non logical container."
    if (.not. allocated(self%l_storage)) error stop &
      "MAC: Error #10: the logical container is not allocated."
    self%l_storage = val
  end subroutine set_gl

  pure subroutine get_al(self, val, at)
    class(container), intent(inout) :: self
    logical, intent(out) :: val
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
    val = self%l_storage(memory_layout)
  end subroutine get_al

  pure subroutine get_ml(self, val, at)
    class(container), intent(inout) :: self
    logical, intent(out) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 1) error stop &
      "MAC: Error #8: trying to get logical value from non logical container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%l_storage)) error stop &
      "MAC: Error #10: the logical container is not allocated."
    val = self%l_storage(at)
  end subroutine get_ml

  function rd_l(self, op, dim)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        logical, intent(in) :: array(:)
        logical :: op
      end function op
    end interface
    integer, intent(in) :: dim

    type(container) :: rd_l

    integer :: shp(self%rank()), lb(self%rank()), arr(self%rank()), &
               outshp(self%rank() - 1), outlb(self%rank() - 1), outarr(self%rank() - 1), &
               vars(self%rank() - 1)

    integer, allocatable :: perm(:, :)

    logical, allocatable :: tmp(:)
    logical :: res

    integer :: i, j
    logical :: dim_passed

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    if ((dim > self%rank()) .or. (dim < 1)) error stop &
      "MAC: Error #14: dim does not reference a valid dimension label."

    shp = self%shape()
    lb = self%lbounds()
    dim_passed = .false.
    do i = 1, self%rank() - 1
      if (i == dim) dim_passed = .true.
      if (dim_passed) then
        outshp(i) = shp(i + 1)
        outlb(i) = lb(i + 1)
        vars(i) = i + 1
      else
        outshp(i) = shp(i)
        outlb(i) = lb(i)
        vars(i) = i
      endif
    enddo

    perm = self%partial_permutation(vars)

    allocate (tmp(shp(dim)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating tmp. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    call rd_l%construct(self%container_type, outshp, outlb, self%lyt)

    do i = 1, size(perm(:, 1))
      arr = perm(i, :)
      dim_passed = .false.
      do j = 1, self%rank() - 1
        if (j == dim) dim_passed = .true.
        if (dim_passed) then
          outarr(j) = arr(j + 1)
        else
          outarr(j) = arr(j)
        endif
      enddo
      do j = 1, shp(dim)
        call self%get(val=tmp(j), at=arr)
      enddo
      res = op(tmp)
      call rd_l%set(val=res, at=outarr)
    enddo

    deallocate (perm)

  end function rd_l

  function ra_l(self, op)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        logical, intent(in) :: array(:)
        logical :: op
      end function
    end interface

    logical :: ra_l

    integer :: shp(self%rank()), sz

    logical, allocatable :: tmp(:, :), hld(:), res(:)

    integer :: dim, stride

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    shp = self%shape()
    sz = self%size()

    allocate (hld(sz), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    hld = self%l_storage

    select case (self%layout())
    case (0)
      do dim = 1, self%rank()

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    case (1)
      do dim = self%rank(), 1, -1

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    end select

    ra_l = hld(1)

    deallocate (hld)

  end function ra_l

  pure subroutine set_ai(self, val, at)
    class(container), intent(inout) :: self
    integer, intent(in) :: val
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
    self%i_storage(memory_layout) = val
  end subroutine set_ai

  pure subroutine set_mi(self, val, at)
    class(container), intent(inout) :: self
    integer, intent(in) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 2) error stop &
      "MAC: Error #8: trying to set integer value to non integer container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%i_storage)) error stop &
      "MAC: Error #10: the integer container is not allocated."
    self%i_storage(at) = val
  end subroutine set_mi

  pure subroutine set_gi(self, val)
    class(container), intent(inout) :: self
    integer, intent(in) :: val
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 2) error stop &
      "MAC: Error #8: trying to set integer value to non integer container."
    if (.not. allocated(self%i_storage)) error stop &
      "MAC: Error #10: the integer container is not allocated."
    self%i_storage = val
  end subroutine set_gi

  pure subroutine get_ai(self, val, at)
    class(container), intent(inout) :: self
    integer, intent(out) :: val
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
    val = self%i_storage(memory_layout)
  end subroutine get_ai

  pure subroutine get_mi(self, val, at)
    class(container), intent(inout) :: self
    integer, intent(out) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 2) error stop &
      "MAC: Error #8: trying to get integer value from non integer container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%i_storage)) error stop &
      "MAC: Error #10: the integer container is not allocated."
    val = self%i_storage(at)
  end subroutine get_mi

  function rd_i(self, op, dim)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        integer, intent(in) :: array(:)
        integer :: op
      end function op
    end interface
    integer, intent(in) :: dim

    type(container) :: rd_i

    integer :: shp(self%rank()), lb(self%rank()), arr(self%rank()), &
               outshp(self%rank() - 1), outlb(self%rank() - 1), outarr(self%rank() - 1), &
               vars(self%rank() - 1)

    integer, allocatable :: perm(:, :)

    integer, allocatable :: tmp(:)
    integer :: res

    integer :: i, j
    logical :: dim_passed

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    if ((dim > self%rank()) .or. (dim < 1)) error stop &
      "MAC: Error #14: dim does not reference a valid dimension label."

    shp = self%shape()
    lb = self%lbounds()
    dim_passed = .false.
    do i = 1, self%rank() - 1
      if (i == dim) dim_passed = .true.
      if (dim_passed) then
        outshp(i) = shp(i + 1)
        outlb(i) = lb(i + 1)
        vars(i) = i + 1
      else
        outshp(i) = shp(i)
        outlb(i) = lb(i)
        vars(i) = i
      endif
    enddo

    perm = self%partial_permutation(vars)

    allocate (tmp(shp(dim)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating tmp. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    call rd_i%construct(self%container_type, outshp, outlb, self%lyt)

    do i = 1, size(perm(:, 1))
      arr = perm(i, :)
      dim_passed = .false.
      do j = 1, self%rank() - 1
        if (j == dim) dim_passed = .true.
        if (dim_passed) then
          outarr(j) = arr(j + 1)
        else
          outarr(j) = arr(j)
        endif
      enddo
      do j = 1, shp(dim)
        call self%get(val=tmp(j), at=arr)
      enddo
      res = op(tmp)
      call rd_i%set(val=res, at=outarr)
    enddo

    deallocate (perm)

  end function rd_i

  function ra_i(self, op)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        integer, intent(in) :: array(:)
        integer :: op
      end function
    end interface

    integer :: ra_i

    integer :: shp(self%rank()), sz

    integer, allocatable :: tmp(:, :), hld(:), res(:)

    integer :: dim, stride

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    shp = self%shape()
    sz = self%size()

    allocate (hld(sz), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    hld = self%i_storage

    select case (self%layout())
    case (0)
      do dim = 1, self%rank()

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    case (1)
      do dim = self%rank(), 1, -1

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    end select

    ra_i = hld(1)

    deallocate (hld)

  end function ra_i

  pure subroutine set_ar(self, val, at)
    class(container), intent(inout) :: self
    real, intent(in) :: val
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
    self%r_storage(memory_layout) = val
  end subroutine set_ar

  pure subroutine set_mr(self, val, at)
    class(container), intent(inout) :: self
    real, intent(in) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 3) error stop &
      "MAC: Error #8: trying to set real value to non real container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%r_storage)) error stop &
      "MAC: Error #10: the real container is not allocated."
    self%r_storage(at) = val
  end subroutine set_mr

  pure subroutine set_gr(self, val)
    class(container), intent(inout) :: self
    real, intent(in) :: val
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 3) error stop &
      "MAC: Error #8: trying to set real value to non real container."
    if (.not. allocated(self%r_storage)) error stop &
      "MAC: Error #10: the real container is not allocated."
    self%r_storage = val
  end subroutine set_gr

  pure subroutine get_ar(self, val, at)
    class(container), intent(inout) :: self
    real, intent(out) :: val
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
    val = self%r_storage(memory_layout)
  end subroutine get_ar

  pure subroutine get_mr(self, val, at)
    class(container), intent(inout) :: self
    real, intent(out) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 3) error stop &
      "MAC: Error #8: trying to get real value from non real container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%r_storage)) error stop &
      "MAC: Error #10: the real container is not allocated."
    val = self%r_storage(at)
  end subroutine get_mr

  function rd_r(self, op, dim)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        real, intent(in) :: array(:)
        real :: op
      end function op
    end interface
    integer, intent(in) :: dim

    type(container) :: rd_r

    integer :: shp(self%rank()), lb(self%rank()), arr(self%rank()), &
               outshp(self%rank() - 1), outlb(self%rank() - 1), outarr(self%rank() - 1), &
               vars(self%rank() - 1)

    integer, allocatable :: perm(:, :)

    real, allocatable :: tmp(:)
    real :: res

    integer :: i, j
    logical :: dim_passed

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    if ((dim > self%rank()) .or. (dim < 1)) error stop &
      "MAC: Error #14: dim does not reference a valid dimension label."

    shp = self%shape()
    lb = self%lbounds()
    dim_passed = .false.
    do i = 1, self%rank() - 1
      if (i == dim) dim_passed = .true.
      if (dim_passed) then
        outshp(i) = shp(i + 1)
        outlb(i) = lb(i + 1)
        vars(i) = i + 1
      else
        outshp(i) = shp(i)
        outlb(i) = lb(i)
        vars(i) = i
      endif
    enddo

    perm = self%partial_permutation(vars)

    allocate (tmp(shp(dim)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating tmp. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    call rd_r%construct(self%container_type, outshp, outlb, self%lyt)

    do i = 1, size(perm(:, 1))
      arr = perm(i, :)
      dim_passed = .false.
      do j = 1, self%rank() - 1
        if (j == dim) dim_passed = .true.
        if (dim_passed) then
          outarr(j) = arr(j + 1)
        else
          outarr(j) = arr(j)
        endif
      enddo
      do j = 1, shp(dim)
        call self%get(val=tmp(j), at=arr)
      enddo
      res = op(tmp)
      call rd_r%set(val=res, at=outarr)
    enddo

    deallocate (perm)

  end function rd_r

  function ra_r(self, op)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        real, intent(in) :: array(:)
        real :: op
      end function
    end interface

    real :: ra_r

    integer :: shp(self%rank()), sz

    real, allocatable :: tmp(:, :), hld(:), res(:)

    integer :: dim, stride

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    shp = self%shape()
    sz = self%size()

    allocate (hld(sz), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    hld = self%r_storage

    select case (self%layout())
    case (0)
      do dim = 1, self%rank()

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    case (1)
      do dim = self%rank(), 1, -1

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    end select

    ra_r = hld(1)

    deallocate (hld)

  end function ra_r

  pure subroutine set_ac(self, val, at)
    class(container), intent(inout) :: self
    complex, intent(in) :: val
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
    self%c_storage(memory_layout) = val
  end subroutine set_ac

  pure subroutine set_mc(self, val, at)
    class(container), intent(inout) :: self
    complex, intent(in) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 4) error stop &
      "MAC: Error #8: trying to set complex value to non complex container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%c_storage)) error stop &
      "MAC: Error #10: the complex container is not allocated."
    self%c_storage(at) = val
  end subroutine set_mc

  pure subroutine set_gc(self, val)
    class(container), intent(inout) :: self
    complex, intent(in) :: val
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 4) error stop &
      "MAC: Error #8: trying to set complex value to non complex container."
    if (.not. allocated(self%c_storage)) error stop &
      "MAC: Error #10: the complex container is not allocated."
    self%c_storage = val
  end subroutine set_gc

  pure subroutine get_ac(self, val, at)
    class(container), intent(inout) :: self
    complex, intent(out) :: val
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
    val = self%c_storage(memory_layout)
  end subroutine get_ac

  pure subroutine get_mc(self, val, at)
    class(container), intent(inout) :: self
    complex, intent(out) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 4) error stop &
      "MAC: Error #8: trying to get complex value from non complex container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%c_storage)) error stop &
      "MAC: Error #10: the complex container is not allocated."
    val = self%c_storage(at)
  end subroutine get_mc

  function rd_c(self, op, dim)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        complex, intent(in) :: array(:)
        complex :: op
      end function op
    end interface
    integer, intent(in) :: dim

    type(container) :: rd_c

    integer :: shp(self%rank()), lb(self%rank()), arr(self%rank()), &
               outshp(self%rank() - 1), outlb(self%rank() - 1), outarr(self%rank() - 1), &
               vars(self%rank() - 1)

    integer, allocatable :: perm(:, :)

    complex, allocatable :: tmp(:)
    complex :: res

    integer :: i, j
    logical :: dim_passed

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    if ((dim > self%rank()) .or. (dim < 1)) error stop &
      "MAC: Error #14: dim does not reference a valid dimension label."

    shp = self%shape()
    lb = self%lbounds()
    dim_passed = .false.
    do i = 1, self%rank() - 1
      if (i == dim) dim_passed = .true.
      if (dim_passed) then
        outshp(i) = shp(i + 1)
        outlb(i) = lb(i + 1)
        vars(i) = i + 1
      else
        outshp(i) = shp(i)
        outlb(i) = lb(i)
        vars(i) = i
      endif
    enddo

    perm = self%partial_permutation(vars)

    allocate (tmp(shp(dim)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating tmp. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    call rd_c%construct(self%container_type, outshp, outlb, self%lyt)

    do i = 1, size(perm(:, 1))
      arr = perm(i, :)
      dim_passed = .false.
      do j = 1, self%rank() - 1
        if (j == dim) dim_passed = .true.
        if (dim_passed) then
          outarr(j) = arr(j + 1)
        else
          outarr(j) = arr(j)
        endif
      enddo
      do j = 1, shp(dim)
        call self%get(val=tmp(j), at=arr)
      enddo
      res = op(tmp)
      call rd_c%set(val=res, at=outarr)
    enddo

    deallocate (perm)

  end function rd_c

  function ra_c(self, op)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        complex, intent(in) :: array(:)
        complex :: op
      end function
    end interface

    complex :: ra_c

    integer :: shp(self%rank()), sz

    complex, allocatable :: tmp(:, :), hld(:), res(:)

    integer :: dim, stride

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    shp = self%shape()
    sz = self%size()

    allocate (hld(sz), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    hld = self%c_storage

    select case (self%layout())
    case (0)
      do dim = 1, self%rank()

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    case (1)
      do dim = self%rank(), 1, -1

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    end select

    ra_c = hld(1)

    deallocate (hld)

  end function ra_c

  pure subroutine set_ardp(self, val, at)
    class(container), intent(inout) :: self
    real(wp), intent(in) :: val
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
    self%rdp_storage(memory_layout) = val
  end subroutine set_ardp

  pure subroutine set_mrdp(self, val, at)
    class(container), intent(inout) :: self
    real(wp), intent(in) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 5) error stop &
      "MAC: Error #8: trying to set real(wp) value to non real(wp) container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%rdp_storage)) error stop &
      "MAC: Error #10: the real(wp) container is not allocated."
    self%rdp_storage(at) = val
  end subroutine set_mrdp

  pure subroutine set_grdp(self, val)
    class(container), intent(inout) :: self
    real(wp), intent(in) :: val
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 5) error stop &
      "MAC: Error #8: trying to set real(wp) value to non real(wp) container."
    if (.not. allocated(self%rdp_storage)) error stop &
      "MAC: Error #10: the real(wp) container is not allocated."
    self%rdp_storage = val
  end subroutine set_grdp

  pure subroutine get_ardp(self, val, at)
    class(container), intent(inout) :: self
    real(wp), intent(out) :: val
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
    val = self%rdp_storage(memory_layout)
  end subroutine get_ardp

  pure subroutine get_mrdp(self, val, at)
    class(container), intent(inout) :: self
    real(wp), intent(out) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 5) error stop &
      "MAC: Error #8: trying to get real(wp) value from non real(wp) container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%rdp_storage)) error stop &
      "MAC: Error #10: the real(wp) container is not allocated."
    val = self%rdp_storage(at)
  end subroutine get_mrdp

  function rd_rdp(self, op, dim)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        real(wp), intent(in) :: array(:)
        real(wp) :: op
      end function op
    end interface
    integer, intent(in) :: dim

    type(container) :: rd_rdp

    integer :: shp(self%rank()), lb(self%rank()), arr(self%rank()), &
               outshp(self%rank() - 1), outlb(self%rank() - 1), outarr(self%rank() - 1), &
               vars(self%rank() - 1)

    integer, allocatable :: perm(:, :)

    real(wp), allocatable :: tmp(:)
    real(wp) :: res

    integer :: i, j
    logical :: dim_passed

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    if ((dim > self%rank()) .or. (dim < 1)) error stop &
      "MAC: Error #14: dim does not reference a valid dimension label."

    shp = self%shape()
    lb = self%lbounds()
    dim_passed = .false.
    do i = 1, self%rank() - 1
      if (i == dim) dim_passed = .true.
      if (dim_passed) then
        outshp(i) = shp(i + 1)
        outlb(i) = lb(i + 1)
        vars(i) = i + 1
      else
        outshp(i) = shp(i)
        outlb(i) = lb(i)
        vars(i) = i
      endif
    enddo

    perm = self%partial_permutation(vars)

    allocate (tmp(shp(dim)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating tmp. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    call rd_rdp%construct(self%container_type, outshp, outlb, self%lyt)

    do i = 1, size(perm(:, 1))
      arr = perm(i, :)
      dim_passed = .false.
      do j = 1, self%rank() - 1
        if (j == dim) dim_passed = .true.
        if (dim_passed) then
          outarr(j) = arr(j + 1)
        else
          outarr(j) = arr(j)
        endif
      enddo
      do j = 1, shp(dim)
        call self%get(val=tmp(j), at=arr)
      enddo
      res = op(tmp)
      call rd_rdp%set(val=res, at=outarr)
    enddo

    deallocate (perm)

  end function rd_rdp

  function ra_rdp(self, op)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        real(wp), intent(in) :: array(:)
        real(wp) :: op
      end function
    end interface

    real(wp) :: ra_rdp

    integer :: shp(self%rank()), sz

    real(wp), allocatable :: tmp(:, :), hld(:), res(:)

    integer :: dim, stride

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    shp = self%shape()
    sz = self%size()

    allocate (hld(sz), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    hld = self%rdp_storage

    select case (self%layout())
    case (0)
      do dim = 1, self%rank()

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    case (1)
      do dim = self%rank(), 1, -1

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    end select

    ra_rdp = hld(1)

    deallocate (hld)

  end function ra_rdp

  pure subroutine set_acdp(self, val, at)
    class(container), intent(inout) :: self
    complex(wp), intent(in) :: val
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
    self%cdp_storage(memory_layout) = val
  end subroutine set_acdp

  pure subroutine set_mcdp(self, val, at)
    class(container), intent(inout) :: self
    complex(wp), intent(in) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 6) error stop &
      "MAC: Error #8: trying to set complex(wp) value to non complex(wp) container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%cdp_storage)) error stop &
      "MAC: Error #10: the complex(wp) container is not allocated."
    self%cdp_storage(at) = val
  end subroutine set_mcdp

  pure subroutine set_gcdp(self, val)
    class(container), intent(inout) :: self
    complex(wp), intent(in) :: val
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 6) error stop &
      "MAC: Error #8: trying to set complex(wp) value to non complex(wp) container."
    if (.not. allocated(self%cdp_storage)) error stop &
      "MAC: Error #10: the complex(wp) container is not allocated."
    self%cdp_storage = val
  end subroutine set_gcdp

  pure subroutine get_acdp(self, val, at)
    class(container), intent(inout) :: self
    complex(wp), intent(out) :: val
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
    val = self%cdp_storage(memory_layout)
  end subroutine get_acdp

  pure subroutine get_mcdp(self, val, at)
    class(container), intent(inout) :: self
    complex(wp), intent(out) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= 6) error stop &
      "MAC: Error #8: trying to get complex(wp) value from non complex(wp) container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%cdp_storage)) error stop &
      "MAC: Error #10: the complex(wp) container is not allocated."
    val = self%cdp_storage(at)
  end subroutine get_mcdp

  function rd_cdp(self, op, dim)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        complex(wp), intent(in) :: array(:)
        complex(wp) :: op
      end function op
    end interface
    integer, intent(in) :: dim

    type(container) :: rd_cdp

    integer :: shp(self%rank()), lb(self%rank()), arr(self%rank()), &
               outshp(self%rank() - 1), outlb(self%rank() - 1), outarr(self%rank() - 1), &
               vars(self%rank() - 1)

    integer, allocatable :: perm(:, :)

    complex(wp), allocatable :: tmp(:)
    complex(wp) :: res

    integer :: i, j
    logical :: dim_passed

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    if ((dim > self%rank()) .or. (dim < 1)) error stop &
      "MAC: Error #14: dim does not reference a valid dimension label."

    shp = self%shape()
    lb = self%lbounds()
    dim_passed = .false.
    do i = 1, self%rank() - 1
      if (i == dim) dim_passed = .true.
      if (dim_passed) then
        outshp(i) = shp(i + 1)
        outlb(i) = lb(i + 1)
        vars(i) = i + 1
      else
        outshp(i) = shp(i)
        outlb(i) = lb(i)
        vars(i) = i
      endif
    enddo

    perm = self%partial_permutation(vars)

    allocate (tmp(shp(dim)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating tmp. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    call rd_cdp%construct(self%container_type, outshp, outlb, self%lyt)

    do i = 1, size(perm(:, 1))
      arr = perm(i, :)
      dim_passed = .false.
      do j = 1, self%rank() - 1
        if (j == dim) dim_passed = .true.
        if (dim_passed) then
          outarr(j) = arr(j + 1)
        else
          outarr(j) = arr(j)
        endif
      enddo
      do j = 1, shp(dim)
        call self%get(val=tmp(j), at=arr)
      enddo
      res = op(tmp)
      call rd_cdp%set(val=res, at=outarr)
    enddo

    deallocate (perm)

  end function rd_cdp

  function ra_cdp(self, op)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        complex(wp), intent(in) :: array(:)
        complex(wp) :: op
      end function
    end interface

    complex(wp) :: ra_cdp

    integer :: shp(self%rank()), sz

    complex(wp), allocatable :: tmp(:, :), hld(:), res(:)

    integer :: dim, stride

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."

    shp = self%shape()
    sz = self%size()

    allocate (hld(sz), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    hld = self%cdp_storage

    select case (self%layout())
    case (0)
      do dim = 1, self%rank()

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    case (1)
      do dim = self%rank(), 1, -1

        allocate (tmp(shp(dim), sz/shp(dim)), res(sz/shp(dim)), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating tmp, res. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        do stride = 1, sz/shp(dim)
          tmp(:, stride) = hld(1 + shp(dim)*(stride - 1):shp(dim)*stride)
          res(stride) = op(tmp(:, stride))
        enddo

        sz = sz/shp(dim)
        deallocate (hld)
        allocate (hld(sz), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "MAC: Error #3: failure allocating hld. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif
        hld = res
        deallocate (tmp, res)

      enddo
    end select

    ra_cdp = hld(1)

    deallocate (hld)

  end function ra_cdp

end module MAC
