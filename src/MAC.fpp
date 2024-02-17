#:set types = ['logical', 'integer', 'real', 'complex', 'real(wp)', 'complex(wp)']
#:set numcode = [1, 2, 3, 4, 5, 6]
#:set storage_name = ['l_storage', 'i_storage', 'r_storage', 'c_storage', 'rdp_storage', 'cdp_storage']
#:set suffixes = ['l', 'i', 'r', 'c', 'rdp', 'cdp']
#:set lst = list(zip(types, numcode, storage_name, suffixes))
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

#:for type, num, storage_name, suffix in lst
  pure subroutine set_a${suffix}$(self, val, at)
    class(container), intent(inout) :: self
    ${type}$, intent(in) :: val
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
    self%${storage_name}$(memory_layout) = val
  end subroutine set_a${suffix}$

  pure subroutine set_m${suffix}$(self, val, at)
    class(container), intent(inout) :: self
    ${type}$, intent(in) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= ${num}$) error stop &
      "MAC: Error #8: trying to set ${type}$ value to non ${type}$ container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%${storage_name}$)) error stop &
      "MAC: Error #10: the ${type}$ container is not allocated."
    self%${storage_name}$(at) = val
  end subroutine set_m${suffix}$

  pure subroutine set_g${suffix}$(self, val)
    class(container), intent(inout) :: self
    ${type}$, intent(in) :: val
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= ${num}$) error stop &
      "MAC: Error #8: trying to set ${type}$ value to non ${type}$ container."
    if (.not. allocated(self%${storage_name}$)) error stop &
      "MAC: Error #10: the ${type}$ container is not allocated."
    self%${storage_name}$ = val
  end subroutine set_g${suffix}$

  pure subroutine get_a${suffix}$(self, val, at)
    class(container), intent(inout) :: self
    ${type}$, intent(out) :: val
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
    val = self%${storage_name}$(memory_layout)
  end subroutine get_a${suffix}$

  pure subroutine get_m${suffix}$(self, val, at)
    class(container), intent(inout) :: self
    ${type}$, intent(out) :: val
    integer, intent(in) :: at
    if (.not. (self%container_initialized)) error stop &
      "MAC: Error #6: container not initialized."
    if (self%container_type /= ${num}$) error stop &
      "MAC: Error #8: trying to get ${type}$ value from non ${type}$ container."
    if ((at < 1) .or. (at > self%size())) error stop &
      "MAC: Error #9: specified memory layout adress is out of bounds."
    if (.not. allocated(self%${storage_name}$)) error stop &
      "MAC: Error #10: the ${type}$ container is not allocated."
    val = self%${storage_name}$(at)
  end subroutine get_m${suffix}$

  function rd_${suffix}$(self, op, dim)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        ${type}$, intent(in) :: array(:)
        ${type}$ :: op
      end function op
    end interface
    integer, intent(in) :: dim

    type(container) :: rd_${suffix}$

    integer :: shp(self%rank()), lb(self%rank()), arr(self%rank()), &
               outshp(self%rank() - 1), outlb(self%rank() - 1), outarr(self%rank() - 1), &
               vars(self%rank() - 1)

    integer, allocatable :: perm(:, :)

    ${type}$, allocatable :: tmp(:)
    ${type}$ :: res

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

    call rd_${suffix}$%construct(self%container_type, outshp, outlb, self%lyt)

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
      call rd_${suffix}$%set(val=res, at=outarr)
    enddo

    deallocate (perm)

  end function rd_${suffix}$

  function ra_${suffix}$(self, op)

    class(container), intent(inout) :: self
    interface
      function op(array)
        import wp
        ${type}$, intent(in) :: array(:)
        ${type}$ :: op
      end function
    end interface

    ${type}$ :: ra_${suffix}$

    integer :: shp(self%rank()), sz

    ${type}$, allocatable :: tmp(:, :), hld(:), res(:)

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

    hld = self%${storage_name}$

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

    ra_${suffix}$ = hld(1)

    deallocate (hld)

  end function ra_${suffix}$

#:endfor
end module MAC
