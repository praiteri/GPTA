module steinhardt_module
  use, intrinsic :: iso_fortran_env, only: real64
  use spherical_harmonics, only: SphericalHarmonics
  use wigner3j_module, only: wigner3j_single
  implicit none
  real(real64), parameter :: PI = 3.1415926535897932384626433832795028841_real64
  private
  public :: Steinhardt

  type :: Wigner3jCache
    integer :: l
    integer :: m1, m2, m3
    real(real64) :: coeff
  end type Wigner3jCache

  type :: Steinhardt
    private
    integer :: m_lmax
    type(SphericalHarmonics) :: m_harmonics
    type(Wigner3jCache), allocatable :: wigner_coefficients(:)
    integer :: num_coeffs
    real(real64), allocatable :: m_q(:)
    real(real64), allocatable :: m_w(:)
    complex(real64), allocatable :: m_qlm(:), m_ylm(:)
  contains
    procedure :: init => Steinhardt_init
    procedure :: compute_q => Steinhardt_compute_q
    procedure :: compute_w => Steinhardt_compute_w
    procedure :: compute_qlm => Steinhardt_compute_qlm
    procedure :: precompute_wigner3j_coefficients => Steinhardt_precompute_wigner3j_coefficients
    procedure :: compute_q_neighbor_products => Steinhardt_compute_q_neighbor_products
    procedure :: size => Steinhardt_size
    procedure :: nlm => Steinhardt_nlm
  end type Steinhardt

contains

  subroutine Steinhardt_init(self, lmax)
  class(Steinhardt), intent(out) :: self
    integer, intent(in) :: lmax
    self%m_lmax = lmax
    call self%m_harmonics%init(lmax)
  end subroutine Steinhardt_init

  function Steinhardt_compute_qlm(self, positions) result(qlm)
  class(Steinhardt), intent(inout) :: self
    real(real64), intent(in) :: positions(:,:)
    complex(real64), allocatable :: qlm(:)
    integer :: num_particles, idx, l, m, i
    complex(real64), allocatable :: ylm(:)
    real(real64) :: pos(3)
    num_particles = size(positions, 2)
    allocate(ylm(self%m_harmonics%nlm()))
    allocate(qlm(self%m_harmonics%nlm()))

    qlm = (0.0_real64, 0.0_real64)
    if (num_particles < 2) return
    do i = 1, num_particles
      pos = positions(:,i)
      pos(:) = pos(:) / norm2(pos)
      call self%m_harmonics%evaluate_pos(pos, ylm)
      idx = 1
      do l = 0, self%m_lmax
        do m = -l, l
        qlm(idx) = qlm(idx) + ylm(idx)
        idx = idx + 1
        end do
      end do
    end do
    qlm = qlm / num_particles
  end function Steinhardt_compute_qlm

  function Steinhardt_compute_q(self, positions) result(q)
  class(Steinhardt), intent(inout) :: self
    real(real64), intent(in) :: positions(:,:)
    real(real64), allocatable :: q(:)
    complex(real64), allocatable :: qlm(:)
    integer :: idx, l, m
    real(real64) :: ql
    allocate(q(self%m_lmax + 1))
    qlm = self%compute_qlm(positions)
    idx = 1
    do l = 0, self%m_lmax
      q(l+1) = 0.0_real64
      ql = 0.0_real64
      do m = -l, l
        ql = ql + abs(qlm(idx) ** 2)
        idx = idx + 1
      end do

      ql = sqrt(4 * PI / (2.0_real64 * l + 1.0_real64) * ql);
      q(l+1) = q(l+1) + ql
    end do
  end function Steinhardt_compute_q

  subroutine Steinhardt_precompute_wigner3j_coefficients(self)
    class(Steinhardt), intent(inout) :: self
    integer :: l, m1, m2, m3
    real(real64) :: w3j

    self%num_coeffs = 0
    do l = 0, self%m_lmax
      do m1 = -l, l
        do m2 = -l, l
          m3 = -m1 - m2
          if (m3 < -l .or. m3 > l) cycle
          w3j = wigner3j_single(real(l,real64), real(l,real64), real(l,real64), &
                                real(m1,real64), real(m2,real64), real(m3,real64))
          if (abs(w3j) > 1e-18_real64) then
            self%num_coeffs = self%num_coeffs + 1
          end if
        end do
      end do
    end do

    allocate(self%wigner_coefficients(self%num_coeffs))
    self%num_coeffs = 0
    do l = 0, self%m_lmax
      do m1 = -l, l
        do m2 = -l, l
          m3 = -m1 - m2
          if (m3 < -l .or. m3 > l) cycle
          w3j = wigner3j_single(real(l,real64), real(l,real64), real(l,real64), &
                                real(m1,real64), real(m2,real64), real(m3,real64))
          if (abs(w3j) > 1e-18_real64) then
            self%num_coeffs = self%num_coeffs + 1
            self%wigner_coefficients(self%num_coeffs) = Wigner3jCache(l, m1, m2, m3, w3j)
          end if
        end do
      end do
    end do
  end subroutine Steinhardt_precompute_wigner3j_coefficients

  function Steinhardt_compute_w(self, positions) result(w)
  class(Steinhardt), intent(inout) :: self
    real(real64), intent(inout) :: positions(:,:)
    real(real64), allocatable :: w(:)
    complex(real64), allocatable :: qlm(:)
    integer :: l, m1, m2, m3, idx1, idx2, idx3, i
    real(real64) :: w3j

    if (.not. allocated(self%wigner_coefficients)) then
      call self%precompute_wigner3j_coefficients
    end if

    allocate(w(self%m_lmax + 1))
    w = 0.0_real64

    qlm = self%compute_qlm(positions)
    do i = 1, self%num_coeffs
      l = self%wigner_coefficients(i)%l
      m1 = self%wigner_coefficients(i)%m1
      m2 = self%wigner_coefficients(i)%m2
      m3 = self%wigner_coefficients(i)%m3
      w3j = self%wigner_coefficients(i)%coeff
      idx1 = l * l + l + m1 + 1
      idx2 = l * l + l + m2 + 1
      idx3 = l * l + l + m3 + 1
      w(l+1) = w(l+1) + real(qlm(idx1) * qlm(idx2) * qlm(idx3) * w3j, real64)
    end do
  end function Steinhardt_compute_w

  function Steinhardt_compute_q_neighbor_products(self, positions, nneigh, lneigh) result(q)
  use moduleDistances  
  class(Steinhardt), intent(inout) :: self
      real(real64), intent(in) :: positions(:,:)
      integer, intent(in) :: nneigh(:)
      integer, intent(in) :: lneigh(:,:)
      real(real64), allocatable :: q(:, :)

      integer :: num_atoms, max_neighbours, i, j, k, idx, l, m
      complex(real64) :: tmp
      complex(real64), allocatable :: qlm(:,:)
      real(real64), allocatable :: rvecs(:, :)
      real(real64) :: norm_i, norm_k, dij(3), dist

      num_atoms = size(positions, 2)
      max_neighbours = size(lneigh, 1)

      if (allocated(qlm)) then
        deallocate(qlm,q,rvecs)
      end if

      allocate(rvecs(3, max_neighbours))
      allocate(q(self%m_lmax + 1, num_atoms))
      q = 0.0_real64
      allocate(qlm(self%m_harmonics%nlm(), num_atoms))

      do i = 1, num_atoms
          if (nneigh(i)==0) cycle
          do j = 1, nneigh(i)
              k = lneigh(j, i)
              dij(1:3) = positions(:, k) - positions(:, i)
              dist = computeDistanceSquaredPBC(dij)
              rvecs(:, j) = dij(:) ! FIX THIS
          end do
          qlm(:, i) = self%compute_qlm(rvecs(:, :nneigh(i)))

      end do
      do i = 1, num_atoms
          if (nneigh(i)==0) cycle
          do j = 1, nneigh(i)
              k = lneigh(j, i)
              idx = 1
              do l = 0, self%m_lmax
                  norm_i = 0.0d0
                  norm_k = 0.0d0
                  tmp = 0.0d0
                  do m = -l, l
                      tmp = tmp + qlm(idx, i) * conjg(qlm(idx, k))
                      norm_i = norm_i + qlm(idx, i) * conjg(qlm(idx, i))
                      norm_k = norm_k + qlm(idx, k) * conjg(qlm(idx, k))
                      idx = idx + 1
                  end do
                  q(l + 1, i) = q(l + 1, i) + real(tmp) / (sqrt(norm_i) * sqrt(norm_k))
              end do
          end do
          q(:, i) = q(:, i) / (nneigh(i))
      end do
      deallocate(rvecs)
      deallocate(qlm)

  end function Steinhardt_compute_q_neighbor_products

  integer function Steinhardt_size(self)
  class(Steinhardt), intent(in) :: self
    Steinhardt_size = self%m_lmax + 1
  end function Steinhardt_size

  integer function Steinhardt_nlm(self)
  class(Steinhardt), intent(in) :: self
    Steinhardt_nlm = self%m_harmonics%nlm()
  end function Steinhardt_nlm

end module steinhardt_module
