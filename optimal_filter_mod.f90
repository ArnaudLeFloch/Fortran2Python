module optimal_filter_mod
  implicit none
  contains

  function optimalfilter3(p_main, p_ref, L) result(p_filtered)
    implicit none
    integer, intent(in) :: L
    real(8), dimension(:), intent(in) :: p_main, p_ref
    real(8), allocatable :: p_filtered(:)
    integer :: Ne, i, j, n
    real(8), allocatable :: g(:), r(:), p_noise(:), f(:)
    real(8), allocatable :: RR(:,:)

    Ne = size(p_main)
    allocate(g(L))
    allocate(r(L))
    allocate(RR(L, L))
    allocate(f(L))
    allocate(p_noise(Ne))
    allocate(p_filtered(Ne))

    ! Compute cross-correlation vector g
    g = 0.0
    do i = 1, L
      do j = 1, Ne
        if (j + i - 1 <= Ne) then
          g(i) = g(i) + p_main(j) * p_ref(j + i - 1)
        end if
      end do
    end do

    ! Compute autocorrelation vector r
    r = 0.0
    do i = 1, L
      do j = 1, Ne
        if (j + i - 1 <= Ne) then
          r(i) = r(i) + p_ref(j) * p_ref(j + i - 1)
        end if
      end do
    end do

    ! Build Toeplitz matrix RR from r
    do i = 1, L
      do j = 1, L
        RR(i, j) = r(abs(i - j) + 1)
      end do
    end do

    ! Solve RR * f = g to get filter coefficients
    call solve_linear_system(RR, g, f, L)

    ! FIR filter - estimate the noise
    p_noise = 0.0
    do n = L + 1, Ne
      do i = 1, L
        p_noise(n) = p_noise(n) + f(i) * p_ref(n - i)
      end do
    end do

    ! Subtract estimated noise from main signal
    do n = 1, Ne
      p_filtered(n) = p_main(n) - p_noise(n)
    end do

    ! Zero out the first L samples
    do n = 1, L
      p_filtered(n) = 0.0
    end do

  end function optimalfilter3

  ! Utility to solve linear system A*x = b using Gauss-Jordan
  subroutine solve_linear_system(A, b, x, N)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N,N), intent(in) :: A
    real(8), dimension(N), intent(in) :: b
    real(8), dimension(N), intent(out) :: x

    real(8), dimension(N,N) :: A_copy
    real(8), dimension(N) :: b_copy
    integer :: ipiv(N), info

    A_copy = A
    b_copy = b

    call dgesv(N, 1, A_copy, N, ipiv, b_copy, N, info)
    if (info /= 0) then
      print *, "Error solving linear system: INFO = ", info
    end if

    x = b_copy

  end subroutine solve_linear_system

end module optimal_filter_mod





