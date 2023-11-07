module io
use global, only: dp, int64
implicit none
private
public readin, readsys, str_padding
contains

! read config
subroutine readin(seed, L0, temp, pres, walkers, burn, prod, snaps, inimode)
    integer, intent(out) :: walkers, burn, prod, snaps
    integer(int64), intent(out) :: seed
    real(dp), intent(out) :: L0, temp, pres
    character(len=*), intent(out) :: inimode
    integer :: iostat
    
    ! set up namelists
    namelist /input/ seed, L0, temp, pres, &
                     walkers, burn, prod, snaps, inimode
    seed = 123456789
    L0 = 10.0
    temp = 1.0
    pres = 1.0
    walkers = 1
    burn = 100
    prod = 100
    snaps = 30
    inimode = 'random'

    ! read
    open(file='input.in', unit=156, iostat=iostat)
    if (iostat .ne. 0) then
    write(*,*) 'Warning: input file opening may have failed.'
    end if

    read(nml=input, unit=156, iostat=iostat)
    if (iostat .ne. 0) then
    write(*,'(a60)') str_padding('Warning: input file reading may have failed. Please', 59)
    write(*,'(a60)') str_padding('double check with the values printed below.', 59)
    end if

    close(156)

    end subroutine readin

! read config
subroutine readsys(nPart, mass, e0, sigma, frcut, asp)
    integer, intent(out) :: nPart
    real(dp), intent(out) :: mass, e0, sigma, frcut
    character(len=*), intent(out) :: asp
    integer :: iostat
    
    ! set up namelists
    namelist /system/  nPart, mass, e0, sigma, frcut, asp
    nPart = 10
    mass = 1.0 
    e0 = 1.0
    sigma = 1.0
    frcut = 0.5
    asp = 'Ar'

    ! read
    open(file='input.in', unit=156, iostat=iostat)
    if (iostat .ne. 0) then
    write(*,*) 'Warning: input file opening may have failed.'
    end if

    read(nml=system, unit=156, iostat=iostat)
    if (iostat .ne. 0) then
    write(*,'(a60)') str_padding('Warning: input file reading may have failed. Please', 59)
    write(*,'(a60)') str_padding('double check with the values printed below.', 59)
    end if

    close(156)

    end subroutine readsys


function str_padding(str, n) result(padded_str)
    integer, intent(in) :: n
    character(len=*), intent(in) :: str
    character(len=n) :: padded_str

    padded_str = trim(str)

    end function str_padding

end module io