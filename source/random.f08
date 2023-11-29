module random
use global, only: dp, int64, pi
implicit none
private
public init_PRNG, randint, randn_matrix
contains


! generate a normal
! using box muller
function randn() result(x)
    implicit none
    real(dp) :: u(2), x

    call random_number(u)
    x = (-2.0_dp*log(u(1)))**0.5_dp*sin(2.0_dp*pi*u(2))
    end function randn

! generate a matrix of normals
! using box muller
function randn_matrix(m, n) result(mat)
    implicit none
    integer, intent(in) :: m, n
    integer :: i, j
    real(dp) :: mat(m,n)

    do i=1,m
    do j=1,n
        mat(i,j) = randn()
    end do
    end do

    end function randn_matrix

! generate a random integer
! from 0 to a (exclusive)
function randint(a) result(i)
    implicit none
    integer, intent(in) :: a
    integer :: i
    real(dp) :: x

    call random_number(x)
    i = a*x

    end function randint

! generate a random binary integer
! valued -1 or +1 with 50/50 chance
function randpm() result(i)
    implicit none
    integer :: i
    real(dp) :: x

    call random_number(x)
    if (x .ge. 0.5_dp) then
        i = 1
        return
    else
        i = -1
        return
    end if

    end function randpm

! initialize PRNG
subroutine init_PRNG(seed)
    implicit none
    integer(int64), intent(in) :: seed
    integer, allocatable :: arr(:)
    integer :: i, n, burn=10000000
    real(dp) :: x
    
    ! get generator size
    call random_seed(size=n)
    allocate(arr(n))

    ! generate more seeds from seed
    do i=1,n
        arr(i) = lcg(seed)
    end do

    ! put array of seeds
    call random_seed(put=arr)

    ! burn in some numbers
    do i=1,burn
        call random_number(x)
    end do

    contains
    ! simple generator used to initialize XOR
    function lcg(s)
        integer :: lcg
        integer(int64) :: s
        if (s == 0) then
            s = 104729
        else
            s = mod(s, 4294967296_int64)
        end if
        s = mod(s * 279470273_int64, 4294967291_int64)
        lcg = int(mod(s, int(huge(0), int64)), kind(0))
        end function lcg
    end subroutine init_PRNG

end module random