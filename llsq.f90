program llsq
  use blas95
  use lapack95

  implicit none

  integer :: m, n, k, i, j
  real, dimension(:,:), allocatable :: b, c, b1, b2, l1
  real, dimension(:,:), allocatable :: y2, y
  real, dimension(:),   allocatable :: d, e 
  real, dimension(:),   allocatable :: tau

  read *, m, n, k
  allocate(d(m))
  allocate(e(k))
  allocate(b(m, n))
  allocate(c(k, n))
  allocate(tau(min(n, k)))
  call read_matrix(b)
  read *, d
  call read_matrix(c)
  read *, e

  call gelqf(c, tau)

  print *, "C = LQ:"
  call print_matrix(c)
  print *, "tau:"
  print *, tau
  
  print *, "c::"
  call print_matrix(c(:, 1:k))
  print *, "e", e

  
  call trsv(c(:, 1:k), e, 'L')
  print *, e
  
  call ormlq(c, tau, b, 'R', 'T')
  print *, "b:"
  call print_matrix(b)

  call gemv(b(:, 1:k), e, d, -1.0)
  print *, "bq1 y1 = d", d
  
  allocate(y2(max(m,k),1))
  y2(:, 1) = d
  print *, y2
  call print_matrix(y2)

  call gels(b(:, k+1:n), y2)

  allocate(y(n, 1))
  y(1:k, 1) = e
  y(k+1:n, 1) = y2(:, 1)

  print *, "y: ", y

  call ormlq(c, tau, y, 'L', 'T')
  
  print *, "wynik: ", y
  ! call print_matrix(b)
  ! print *, d

  
contains
  subroutine read_matrix(a)
    integer :: i, j
    real, dimension(:,:) :: a

    do i = 1, size(a, 1)
       read *, (a(i, j), j = 1, size(a, 2))
    end do
  end subroutine read_matrix

  subroutine print_matrix(a)
    integer :: i, j
    real, dimension(:,:) :: a

    do i = 1, size(a, 1)
       print *, (a(i, j), j = 1, size(a, 2))
    end do
  end subroutine print_matrix
end program llsq
