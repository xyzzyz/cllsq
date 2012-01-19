program llsq
  use blas95
  use lapack95

  implicit none

  integer :: n
  real, dimension(:,:), allocatable :: b, c
  real, dimension(:,:), allocatable :: y
  real, dimension(:),   allocatable :: d, e 

  call read_cllsq(b, d, c, e)
  
  n = size(b, 2)
  allocate(y(n, 1))

  call solve_cllsq(b, d, c, e, y)

  print *, "wynik: ", y

contains
  subroutine read_cllsq(b, d, c, e) 
    integer :: m, n, k
    real, dimension(:,:), allocatable, intent(out) :: b, c
    real, dimension(:),   allocatable, intent(out) :: d, e 

    read *, m, k, n

    allocate(d(m))
    allocate(e(k))
    allocate(b(m, n))
    allocate(c(k, n))
    call read_matrix(b)
    call read_matrix(c)
    read *, d
    read *, e

  end subroutine read_cllsq

  subroutine solve_cllsq(b, d, c, e, y)
    real, dimension(:,:), intent(inout) :: b, c
    real, dimension(:,:), intent(inout) :: y
    real, dimension(:),   intent(inout) :: d, e 
    real, dimension(:,:), allocatable :: y2
    real, dimension(:),   allocatable :: tau
    integer :: m, n, k
    
    m = size(d)
    k = size(e)
    n = size(b, 2)

    allocate(tau(min(n, k)))


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

    y(1:k, 1) = e
    y(k+1:n, 1) = y2(1:n-k, 1)
    
    ! oblicz x = Q y
    call ormlq(c, tau, y, 'L', 'T')

  end subroutine solve_cllsq
  subroutine read_matrix(a)
    integer :: i, j
    real, dimension(:,:), intent(out) :: a

    do i = 1, size(a, 1)
       read *, (a(i, j), j = 1, size(a, 2))
    end do
  end subroutine read_matrix

  subroutine print_matrix(a)
    integer :: i, j
    real, dimension(:,:), intent(in) :: a

    do i = 1, size(a, 1)
       print *, (a(i, j), j = 1, size(a, 2))
    end do
  end subroutine print_matrix
end program llsq
