program main
  use optimal_filter_mod
  implicit none
  integer :: Ne, L, i
  real(8), allocatable :: p_main(:), p_ref(:), p_filtered(:)

  ! Signal length and filter order
  Ne = 1000
  L = 50

  ! Allocate and generate example signals
  allocate(p_main(Ne))
  allocate(p_ref(Ne))

  ! Create a reference signal (e.g., a sine wave)
  do i = 1, Ne
    p_ref(i) = sin(2.0d0 * 3.141592653589d0 * i / 50.0d0)
  end do

  ! Create main signal: reference + some noise
  call random_seed()
  do i = 1, Ne
    call random_number(p_main(i))
    p_main(i) = p_ref(i) + 0.5d0 * (p_main(i) - 0.5d0)  ! Add noise
  end do

  ! Call the optimal filter
  p_filtered = optimalfilter3(p_main, p_ref, L)

  ! Print first few values for comparison
  print *, "Index   p_main    p_filtered"
  do i = 1, 10
    print "(I5, 3X, F8.4, 3X, F8.4)", i, p_main(i), p_filtered(i)
  end do
  
  ! Open a file for writing
  open(unit=10, file='p_main.txt', status='unknown')
  ! Write data
  do i = 1, Ne
    write(10, *) p_main(i)
  end do
  ! Close the file
  close(10)
  
  ! Open a file for writing
  open(unit=10, file='p_filtered.txt', status='unknown')
  ! Write data
  do i = 1, Ne
    write(10, *) p_filtered(i)
  end do
  ! Close the file
  close(10)

  ! Cleanup
  deallocate(p_main, p_ref, p_filtered)

end program main
