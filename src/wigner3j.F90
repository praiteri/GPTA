module wigner3j_module
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: wigner3j, wigner3j_single

  real(real64), parameter :: eps = epsilon(1.0_real64)
  real(real64), parameter :: tiny_value = 1e-20 !tiny(1.0_real64)
  real(real64), parameter :: huge_value = huge(1.0_real64)
  real(real64), parameter :: tolerance = 1.0e-10_real64

contains

  function wigner3j(l2, l3, m1, m2, m3) result(vec)
    real(real64), intent(in) :: l2, l3, m1, m2, m3
    real(real64), allocatable :: vec(:)
    real(real64) :: huge_sentinel, tiny_sentinel
    real(real64) :: l1min, l1max, l1, an, prev_an, beta
    integer :: size, i, j, k, s
    real(real64) :: l1_mid_minus_one, l1_mid, l1_mid_plus_one, lambda
    real(real64) :: tot, c1
    logical :: av

    huge_sentinel = sqrt(sqrt(huge_value / 20.0_real64))
    tiny_sentinel = sqrt(tiny_value * 2)

    if (.not. (abs(m1 + m2 + m3) < eps .and. &
               abs(m2) <= l2 + eps .and. &
               abs(m3) <= l3 + eps)) then
      allocate(vec(1))
      vec = 0.0_real64
      return
    end if

    l1min = max(abs(l2 - l3), abs(m1))
    l1max = l2 + l3

    size = int(floor(l1max - l1min + 1.0_real64 + eps))
    allocate(vec(size))
    vec = 0.0_real64

    if (size == 1) then
      vec(1) = (-1.0_real64)**int(abs(l2 + m2 - l3 + m3)) / sqrt(l1min + l2 + l3 + 1.0_real64)
    else
      vec(1) = tiny_sentinel
      l1 = l1min

      if (abs(l1min) < tolerance) then
        an = alpha_term_l1_zero(l1, l2, l3, m1, m2, m3)
      else
        an = alpha_term(l1min, l2, l3, m1, m2, m3)
      end if

      vec(2) = an * vec(1)

      if (size > 2) then
        vec(1) = tiny_sentinel

        if (abs(l1min) < tolerance) then
          an = alpha_term_l1_zero(l1, l2, l3, m1, m2, m3)
        else
          an = alpha_term(l1min, l2, l3, m1, m2, m3)
        end if

        vec(2) = an * vec(1)

        i = 2
        av = .false.
        do
          i = i + 1
          prev_an = an
          l1 = l1 + 1.0_real64

          an = alpha_term(l1, l2, l3, m1, m2, m3)
          beta = beta_term(l1, l2, l3, m1)

          vec(i) = an * vec(i - 1) + beta * vec(i - 2)

          if (abs(vec(i)) > huge_sentinel) then
            vec(1:i) = vec(1:i) / huge_sentinel
          end if

          if (av) exit
          if (abs(an) - abs(prev_an) > 0.0_real64) av = .true.
          if (i >= size) exit
        end do

        if (i /= size) then
          l1_mid_minus_one = vec(i - 2)
          l1_mid = vec(i - 1)
          l1_mid_plus_one = vec(i)

          vec(size) = tiny_sentinel

          l1 = l1max
          an = alpha_term_backward(l1, l2, l3, m1, m2, m3)
          vec(size - 1) = an * vec(size)

          j = size - 1
          do
            j = j - 1
            l1 = l1 - 1.0_real64

            an = alpha_term_backward(l1, l2, l3, m1, m2, m3)
            beta = beta_term_backward(l1, l2, l3, m1)

            vec(j) = an * vec(j + 2) + beta * vec(j + 2)

            if (abs(vec(j)) > huge_sentinel) then
              vec(j:size) = vec(j:size) / huge_sentinel
            end if

            if (j <= i - 2) exit
          end do

          lambda = (l1_mid_plus_one * vec(j + 2) + l1_mid * vec(j + 1) + l1_mid_minus_one * vec(j)) / &
            (l1_mid_plus_one * l1_mid_plus_one + l1_mid * l1_mid + l1_mid_minus_one * l1_mid_minus_one)

          vec(1:j-1) = vec(1:j-1) * lambda
        end if
      end if
    end if

    tot = 0.0_real64
    do k = 1, size
      tot = tot + (2.0_real64 * (l1min + real(k - 1, real64)) + 1.0_real64) * (vec(k)**2)
    end do

    s = merge(-1, 1, vec(size) < 0.0_real64)
    c1 = (-1.0_real64)**(int(l2 - l3 - m1)) * s
    vec = vec * c1 / sqrt(tot)

  end function wigner3j

  function wigner3j_single(l1, l2, l3, m1, m2, m3) result(val)
    real(real64), intent(in) :: l1, l2, l3, m1, m2, m3
    real(real64) :: val
    real(real64), allocatable :: r(:)
    real(real64) :: l1min
    integer :: index

    if (.not. (abs(m1 + m2 + m3) < tolerance .and. &
               abs(floor(l1 + l2 + l3) - (l1 + l2 + l3)) < tolerance .and. &
               l3 >= abs(l1 - l2) .and. &
               l3 <= l1 + l2 .and. &
               abs(m1) <= l1 .and. &
               abs(m2) <= l2 .and. &
               abs(m3) <= l3)) then
      val = 0.0_real64
      return
    end if

    l1min = max(abs(l2 - l3), abs(m1))
    index = int(l1 - l1min) + 1

    allocate(r(index))  ! Allocate the 'r' array
    r = wigner3j(l2, l3, m1, m2, m3)
    val = r(index)
    deallocate(r)  ! Deallocate the 'r' array

  end function wigner3j_single

  pure function alpha_term_l1_zero(l1, l2, l3, m1, m2, m3) result(val)
    real(real64), intent(in) :: l1, l2, l3, m1, m2, m3
    real(real64) :: val
    val = -(m3 - m2 + 2.0_real64 * w3b(l1, l2, l3, m1, m2, m3)) / w3a(1.0_real64, l2, l3, m1)
  end function alpha_term_l1_zero

  pure function alpha_term(l1, l2, l3, m1, m2, m3) result(val)
    real(real64), intent(in) :: l1, l2, l3, m1, m2, m3
    real(real64) :: val
    val = -w3b(l1, l2, l3, m1, m2, m3) / (l1 * w3a(l1 + 1.0_real64, l2, l3, m1))
  end function alpha_term

  pure function alpha_term_backward(l1, l2, l3, m1, m2, m3) result(val)
    real(real64), intent(in) :: l1, l2, l3, m1, m2, m3
    real(real64) :: val
    val = -w3b(l1, l2, l3, m1, m2, m3) / ((l1 + 1.0_real64) * w3a(l1, l2, l3, m1))
  end function alpha_term_backward

  pure function beta_term(l1, l2, l3, m1) result(val)
    real(real64), intent(in) :: l1, l2, l3, m1
    real(real64) :: val
    val = -(l1 + 1.0_real64) * w3a(l1, l2, l3, m1) / (l1 * w3a(l1 + 1.0_real64, l2, l3, m1))
  end function beta_term

  pure function beta_term_backward(l1, l2, l3, m1) result(val)
    real(real64), intent(in) :: l1, l2, l3, m1
    real(real64) :: val
    val = -l1 * w3a(l1 + 1.0_real64, l2, l3, m1) / ((l1 + 1.0_real64) * w3a(l1, l2, l3, m1))
  end function beta_term_backward

  pure function w3a(l1, l2, l3, m1) result(val)
    real(real64), intent(in) :: l1, l2, l3, m1
    real(real64) :: val
    real(real64) :: t1, t2, t3
    t1 = l1**2 - (l2 - l3)**2
    t2 = (l2 + l3 + 1.0_real64)**2 - l1**2
    t3 = l1**2 - m1**2
    val = sqrt(t1 * t2 * t3)
  end function w3a

  pure function w3b(l1, l2, l3, m1, m2, m3) result(val)
    real(real64), intent(in) :: l1, l2, l3, m1, m2, m3
    real(real64) :: val
    real(real64) :: t1, t2, t3, t4
    t1 = -(2.0_real64 * l1 + 1.0_real64)
    t2 = l2 * (l2 + 1.0_real64) * m1
    t3 = l3 * (l3 + 1.0_real64) * m1
    t4 = l1 * (l1 + 1.0_real64) * (m3 - m2)
    val = t1 * (t2 - t3 - t4)
  end function w3b

end module wigner3j_module
