module legendre
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: AssocLegendreP

  type AssocLegendreP
    private
    integer :: l_max = 0
    real(real64), allocatable :: m_a(:,:)
    real(real64), allocatable :: m_b(:,:)
    real(real64), allocatable :: m_cache(:,:)
  contains
    procedure :: init => AssocLegendreP_init
    procedure :: evaluate => AssocLegendreP_evaluate
    procedure :: evaluate_batch_vec => AssocLegendreP_evaluate_batch_vec
    procedure :: evaluate_batch_sub => AssocLegendreP_evaluate_batch_sub
    procedure :: work_array => AssocLegendreP_work_array
    procedure :: work_array_size => AssocLegendreP_work_array_size
  end type AssocLegendreP

  interface AssocLegendreP
    module procedure AssocLegendreP_constructor
  end interface AssocLegendreP

contains

  type(AssocLegendreP) function AssocLegendreP_constructor(lm) result(self)
    integer, intent(in) :: lm
    call self%init(lm)
  end function AssocLegendreP_constructor

  subroutine AssocLegendreP_init(self, lm)
  class(AssocLegendreP), intent(out) :: self
    integer, intent(in) :: lm
    integer :: m, l
    self%l_max = lm
    allocate(self%m_a(0:lm, 0:lm))
    allocate(self%m_b(0:lm, 0:lm))
    allocate(self%m_cache(0:lm, 0:lm))

    do m = 0, lm
      self%m_a(m, m) = amm(m)
      do l = m + 1, lm
        self%m_a(l, m) = alm(l, m)
        self%m_b(l, m) = blm(l, m)
      end do
    end do
  end subroutine AssocLegendreP_init

  recursive real(real64) function AssocLegendreP_evaluate(self, l, m, x) result(val)
  class(AssocLegendreP), intent(in) :: self
    integer, intent(in) :: l, m
    real(real64), intent(in) :: x
    if (m == l) then
      val = self%m_a(l, m) * (1.0_real64 - x * x) ** (0.5_real64 * m)
    else if (m + 1 == l) then
      val = self%m_a(l, m) * x * self%evaluate(m, m, x)
    else
      val = self%m_a(l, m) * x * self%evaluate(l - 1, m, x) + &
        self%m_b(l, m) * self%evaluate(l - 2, m, x)
    end if
  end function AssocLegendreP_evaluate

  subroutine AssocLegendreP_evaluate_batch_sub(self, x, result)
  class(AssocLegendreP), intent(inout) :: self
    real(real64), intent(in) :: x
    real(real64), intent(out) :: result(:)
    integer :: idx, m, l
    real(real64) :: tmp1, sqrt_tmp1
    idx = 1
    tmp1 = 1.0_real64 - x * x
    sqrt_tmp1 = sqrt(tmp1)

    do m = 0, self%l_max
      do l = m, self%l_max
        if (l == m) then
          result(idx) = self%m_a(l, m) * (1.0_real64 - x * x) ** (0.5_real64 * m)
        else if (l == (m + 1)) then
          result(idx) = self%m_a(l, m) * x * self%m_cache(l - 1, m)
        else
          result(idx) = self%m_a(l, m) * x * self%m_cache(l - 1, m) + &
            self%m_b(l, m) * self%m_cache(l - 2, m)
        end if
        self%m_cache(l, m) = result(idx)
        idx = idx + 1
      end do
    end do
  end subroutine AssocLegendreP_evaluate_batch_sub

  function AssocLegendreP_work_array(self) result(arr)
  class(AssocLegendreP), intent(in) :: self
    real(real64), allocatable :: arr(:)
    allocate(arr((self%l_max + 1) * (self%l_max + 2) / 2))
  end function AssocLegendreP_work_array

  integer function AssocLegendreP_work_array_size(self)
  class(AssocLegendreP), intent(in) :: self
    AssocLegendreP_work_array_size = (self%l_max + 1) * (self%l_max + 2) / 2
  end function AssocLegendreP_work_array_size

  function AssocLegendreP_evaluate_batch_vec(self, x) result(arr)
  class(AssocLegendreP), intent(inout) :: self
    real(real64), intent(in) :: x
    real(real64), allocatable :: arr(:)
    arr = self%work_array()
    call self%evaluate_batch_sub(x, arr)
  end function AssocLegendreP_evaluate_batch_vec

  real(real64) function amm(m)
    integer, intent(in) :: m
    real(real64), parameter :: pi4 = 4.0_real64 * acos(-1.0_real64)
    real(real64) :: result
    integer :: k
    result = 1.0_real64
    do k = 1, m
      result = result * (2.0_real64 * k + 1.0_real64) / (2.0_real64 * k)
    end do
    amm = sqrt(result / pi4)
  end function amm

  real(real64) function alm(l, m)
    integer, intent(in) :: l, m
    alm = sqrt((4.0_real64 * l * l - 1.0_real64) / (1.0_real64 * l * l - m * m))
  end function alm

  real(real64) function blm(l, m)
    integer, intent(in) :: l, m
    blm = -sqrt((2.0_real64 * l + 1.0_real64) * ((l - 1.0_real64) * (l - 1.0_real64) - m * m) / &
      ((2.0_real64 * l - 3.0_real64) * (1.0_real64 * l * l - m * m)))
  end function blm

end module legendre
