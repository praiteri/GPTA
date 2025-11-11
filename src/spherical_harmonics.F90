module spherical_harmonics
  use iso_fortran_env, only: real64
  use legendre, only: AssocLegendreP
  implicit none
  private
  public :: SphericalHarmonics

  type :: SphericalHarmonics
    private
    logical :: m_phase = .true.
    integer :: m_lmax = 0
    type(AssocLegendreP) :: m_plm_evaluator
    real(real64), allocatable :: m_plm(:)
  contains
    procedure :: init => SphericalHarmonics_init
    procedure :: evaluate_cvec_theta_phi => SphericalHarmonics_evaluate_cvec_theta_phi
    procedure :: evaluate_cvec_pos => SphericalHarmonics_evaluate_cvec_pos
    procedure :: evaluate_theta_phi => SphericalHarmonics_evaluate_theta_phi
    procedure :: evaluate_pos => SphericalHarmonics_evaluate_pos
    procedure :: nlm => SphericalHarmonics_nlm
  end type SphericalHarmonics

  interface SphericalHarmonics
    module procedure SphericalHarmonics_constructor
  end interface SphericalHarmonics

contains

  type(SphericalHarmonics) function SphericalHarmonics_constructor(lm, phase) result(self)
    integer, intent(in) :: lm
    logical, intent(in), optional :: phase
    if (present(phase)) then
      self%m_phase = phase
    end if
    self%m_lmax = lm
    self%m_plm_evaluator = AssocLegendreP(lm)
    allocate(self%m_plm(self%m_plm_evaluator%work_array_size()))
  end function SphericalHarmonics_constructor

  subroutine SphericalHarmonics_init(self, lm, phase)
  class(SphericalHarmonics), intent(out) :: self
    integer, intent(in) :: lm
    logical, intent(in), optional :: phase
    if (present(phase)) then
      self%m_phase = phase
    end if
    self%m_lmax = lm
    self%m_plm_evaluator = AssocLegendreP(lm)
    allocate(self%m_plm(self%m_plm_evaluator%work_array_size()))
  end subroutine SphericalHarmonics_init

  function SphericalHarmonics_evaluate_cvec_theta_phi(self, theta, phi) result(cvec)
  class(SphericalHarmonics), intent(inout) :: self
    real(real64), intent(in) :: theta, phi
    complex(real64), allocatable :: cvec(:)
    allocate(cvec(self%nlm()))
    call self%evaluate_theta_phi(theta, phi, cvec)
  end function SphericalHarmonics_evaluate_cvec_theta_phi

  subroutine SphericalHarmonics_evaluate_theta_phi(self, theta, phi, result)
  class(SphericalHarmonics), intent(inout) :: self
    real(real64), intent(in) :: theta, phi
    complex(real64), intent(out) :: result(:)
    real(real64) :: ct
    complex(real64) :: c, cm
    integer :: l, m, l_offset, plm_idx, sign
    ct = cos(theta)
    call self%m_plm_evaluator%evaluate_batch_sub(ct, self%m_plm)
    plm_idx = 1
    do l = 0, self%m_lmax
      l_offset = l * (l + 1) + 1
      result(l_offset) = self%m_plm(plm_idx)
      plm_idx = plm_idx + 1
    end do
    c = exp(cmplx(0.0_real64, phi, real64))
    cm = c
    do m = 1, self%m_lmax
      sign = merge(-1, 1, self%m_phase .and. (mod(m, 2) == 1))
      do l = m, self%m_lmax
        l_offset = l * (l + 1) + 1
        result(l_offset - m) = sign * self%m_plm(plm_idx) * conjg(cm)
        result(l_offset + m) = sign * self%m_plm(plm_idx) * cm
        if (mod(m, 2) == 1) result(l_offset - m) = -result(l_offset - m)
        plm_idx = plm_idx + 1
      end do
      cm = cm * c
    end do
  end subroutine SphericalHarmonics_evaluate_theta_phi

  function SphericalHarmonics_evaluate_cvec_pos(self, pos) result(cvec)
  class(SphericalHarmonics), intent(inout) :: self
    real(real64), intent(in) :: pos(3)
    complex(real64), allocatable :: cvec(:)
    allocate(cvec(self%nlm()))
    call self%evaluate_pos(pos, cvec)
  end function SphericalHarmonics_evaluate_cvec_pos

  subroutine SphericalHarmonics_evaluate_pos(self, pos, result)
  class(SphericalHarmonics), intent(inout) :: self
    real(real64), intent(in) :: pos(3)
    complex(real64), intent(out) :: result(:)
    real(real64), parameter :: epsilon = 1e-12_real64
    real(real64) :: ct, st
    complex(real64) :: c, cm
    integer :: l, m, l_offset, plm_idx, sign
    ct = pos(3)
    call self%m_plm_evaluator%evaluate_batch_sub(ct, self%m_plm)
    st = merge(sqrt(1.0_real64 - ct * ct), 0.0_real64, abs(1.0_real64 - ct) > epsilon)
    plm_idx = 1
    do l = 0, self%m_lmax
      l_offset = l * (l + 1) + 1
      result(l_offset) = self%m_plm(plm_idx)
      plm_idx = plm_idx + 1
    end do
    c = cmplx(merge(pos(1) / st, 0.0_real64, st > epsilon), &
      merge(pos(2) / st, 0.0_real64, st > epsilon), real64)
    cm = c
    do m = 1, self%m_lmax
      sign = merge(-1, 1, self%m_phase .and. (mod(m, 2) == 1))
      do l = m, self%m_lmax
        l_offset = l * (l + 1) + 1
        result(l_offset - m) = sign * self%m_plm(plm_idx) * conjg(cm)
        result(l_offset + m) = sign * self%m_plm(plm_idx) * cm
        if (mod(m, 2) == 1) result(l_offset - m) = -result(l_offset - m)
        plm_idx = plm_idx + 1
      end do
      cm = cm * c
    end do
  end subroutine SphericalHarmonics_evaluate_pos

  integer function SphericalHarmonics_nlm(self)
  class(SphericalHarmonics), intent(in) :: self
    SphericalHarmonics_nlm = (self%m_lmax + 1) * (self%m_lmax + 1)
  end function SphericalHarmonics_nlm

end module spherical_harmonics
