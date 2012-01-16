program llsq
  use blas95
  use lapack95

  implicit none

  integer :: m, n, k
  real, dimension(:,:), allocatable :: b, c
  real, dimension(:,:), allocatable :: y2, y
  real, dimension(:),   allocatable :: d, e 
  real, dimension(:),   allocatable :: tau

  read *, m, k, n
  allocate(d(m))
  allocate(e(k))
  allocate(b(m, n))
  allocate(c(k, n))
  allocate(tau(min(n, k)))
  call read_matrix(b)
  call read_matrix(c)
  read *, d
  read *, e
  
  ! dokonaj rozkładu C = LQ
  call gelqf(c, tau)
  
  ! rozwiąż układ L_1 y_1 = e
  call trsv(c(:, 1:k), e, 'L')

  ! oblicz BQ
  call ormlq(c, tau, b, 'R', 'T')

  ! oblicz d - B Q_1 

  call gemv(b(:, 1:k), e, d, -1.0, 1.0)

  allocate(y2(max(m,k),1))
  y2(:, 1) = d

  ! rozwiąż RLZNK dla B Q_2 y_2 = d - B Q_1
  call gels(b(:, k+1:n), y2)

  allocate(y(n, 1))
  y(1:k, 1) = e
  y(k+1:n, 1) = y2(1:n-k, 1)

  ! oblicz x = Q y
  call ormlq(c, tau, y, 'L', 'T')
  
  print *, "wynik: ", y
  
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
