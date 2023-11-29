module io
use global, only: dp, int64
implicit none
private
public mc_readin, md_readin, readsys
contains

! read MC config
subroutine mc_readin(seed, L0, temp, press, xvar, vvar, burn, prod, snaps, inimode)
    integer, intent(out) :: burn, prod, snaps
    integer(int64), intent(out) :: seed
    real(dp), intent(out) :: L0, temp, press, xvar, vvar
    character(len=*), intent(out) :: inimode
    integer :: iostat
    
    ! set up namelists
    namelist /input/ seed, L0, temp, press, xvar, vvar, &
                       burn, prod, snaps, inimode
    seed = 123456789
    L0 = 10.0
    temp = 1.0
    press = 1.0
    xvar = 0.005
    vvar = 0.005
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
    write(*,'(a45)') 'Warning: input file reading may have failed.'
    write(*,'(a48)') 'Please double check the values printed below...'
    end if

    close(156)

    end subroutine mc_readin

! read MD config
subroutine md_readin(seed, L0, dt, temp, press, gamma, steps, snaps, inimode)
    integer, intent(out) :: steps, snaps
    integer(int64), intent(out) :: seed
    real(dp), intent(out) :: L0, dt, temp, press, gamma
    character(len=*), intent(out) :: inimode
    integer :: iostat
    
    ! set up namelists
    namelist /input/ seed, L0, dt, temp, press, gamma, steps, snaps, inimode
    seed = 123456789
    L0 = 10.0
    dt = 0.001
    temp = 1.0
    press = 1.0
    gamma = 0.0
    steps = 100
    snaps = 30
    inimode = 'random'

    ! read
    open(file='input.in', unit=156, iostat=iostat)
    if (iostat .ne. 0) then
    write(*,*) 'Warning: input file opening may have failed.'
    end if

    read(nml=input, unit=156, iostat=iostat)
    if (iostat .ne. 0) then
    write(*,'(a45)') 'Warning: input file reading may have failed.'
    write(*,'(a48)') 'Please double check the values printed below...'
    end if

    close(156)

    end subroutine md_readin

! read config
subroutine readsys(npart, mass, e0, sigma, rcut, asp)
    integer, intent(out) :: npart
    real(dp), intent(out) :: mass, e0, sigma, rcut
    character(len=*), intent(out) :: asp
    integer :: iostat
    
    ! set up namelists
    namelist /system/  npart, mass, e0, sigma, rcut, asp
    npart = 10
    mass = 1.0 
    e0 = 1.0
    sigma = 1.0
    rcut = 3.0
    asp = 'Ar'

    ! read
    open(file='input.in', unit=156, iostat=iostat)
    if (iostat .ne. 0) then
    write(*,*) 'Warning: input file opening may have failed.'
    end if

    read(nml=system, unit=156, iostat=iostat)
    if (iostat .ne. 0) then
    write(*,'(a45)') 'Warning: input file reading may have failed.'
    write(*,'(a48)') 'Please double check the values printed below...'
    end if

    close(156)

    end subroutine readsys

end module io