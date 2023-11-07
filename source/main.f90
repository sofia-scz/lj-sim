program main
    use global, only: dp, xp, int64, zero, & 
                      amuang_to_kglt, evang_to_bar, &
                      nPart, mass, e0, sigma, frcut
    use random, only: init_PRNG
    use io, only: readin, readsys, str_padding
    use mcmc, only: mcmc_move_NPT
    use potential, only: compute_poten, compute_vpress, pbc_table
    use omp_lib
    implicit none

    real(dp), allocatable :: x(:,:)
    real(dp) :: L, temp, pres, t0, tf, energy, & 
                    volume, density, vpress
    integer :: i, n, walkers, burn, prod, snaps
    integer(int64) :: seed
    character(len=2) :: dummy, asp
    character(len=10) :: inimode
    character(len=14) :: fstate='finalstate.out', fpos='positions.in'

    ! call init
    call init()

    ! set up variables
    call setup_vars()
    call flush()

    ! start computations
    write(*,*) ''
    write(*,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    write(*,*) '~~~~~~~~~~~~ Starting computations ~~~~~~~~~~~~'
    write(*,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    write(*,*) ''
    call flush()

    t0 = omp_get_wtime()
    write(*,*) 'Starting burn in...'
    call flush()

    do n=1,burn
        call mcmc_move_NPT(nPart, temp, pres, L, x, energy)
        volume = L**3.0_xp
        density = nPart/volume*amuang_to_kglt*mass
        vpress = compute_vpress(nPart, temp, L, x)*evang_to_bar
        call bin_step()

        ! force data dump every 1 per cent
        if (mod(i-mod(i,prod/snaps),prod/snaps).eq.0) then
            call flush()
        end if
    end do

    tf = omp_get_wtime()
    write(*,'(a34,f12.4,a4)') str_padding('Burn in finished! Time spent: '&
                                                        ,33), tf-t0, 'sec'

    t0 = omp_get_wtime()
    write(*,*) ''
    write(*,*) 'Starting production...'

    open(file='table.out', unit=236)
    write(236,'(4a18)') 'Energy', 'BoxLength', 'Density', 'VirialPress'
    open(unit=145, file='snapshots.out')

    do n=1,prod
        call mcmc_move_NPT(nPart, temp, pres, L, x, energy)
        volume = L**3.0_xp
        density = nPart/volume*amuang_to_kglt
        vpress = compute_vpress(nPart, temp, L, x)*evang_to_bar
        write(236,'(e18.8,2f18.8,e18.8)') energy/nPart, L, density, vpress

        ! write snapshots
        if (mod(n-mod(prod,int(prod/snaps)),int(prod/snaps)).eq.0) then
            write(145,*) nPart
            write(145,'(a12, 9f6.2, a2)') 'Lattice= "', L, .0, .0, .0, L, .0, .0, .0, L, '"'
            do i=1,nPart
                write(145,'(a4, 3f16.8)') asp, x(i,:)
            end do
            call flush()
        else if (mod(n-mod(prod,int(prod/100)),int(prod/100)).eq.0) then ! force data dump every 1 per cent
            call flush()
        end if
    end do

    close(145)
    close(236)
    tf = omp_get_wtime()

    write(*,'(a34,f12.4,a4)') str_padding('Production finished! Time spent: '&
                                        ,33), tf-t0, 'sec'
    write(*,*) ''
    write(*,*) 'Saving final state data...'

    call save_fstep()

    write(*,*) 'Success!'
    write(*,*) 'Execution complete. '

contains

! initialize program
subroutine init()
    implicit none
    ! start writing output
    write(*,*) ''
    write(*,*) '-------------------------------------------------------------'
    write(*,*) '| Welcome to my LJ fluid MCMC simulation!!                  |'
    write(*,*) '|                                                           |'
    write(*,*) '| This program was written by C. Dacal and S. Scozziero for |'
    write(*,*) '| the course "Introduccion a la Simulacion Computacional"   |'
    write(*,*) '| taught by Dr. Pastorino at Instituto Sabato.              |'
    write(*,*) '|                                                           |'
    write(*,*) '| We hope you enjoy using our program (✿◠‿◠)                |'
    write(*,*) '-------------------------------------------------------------'

    write(*,*) ''
    write(*,*) 'Reading input files...'
    call readin(seed, L, temp, pres, walkers, burn, prod, snaps, inimode)
    call readsys(nPart, mass, e0, sigma, frcut, asp)
    write(*,*) 'Success!'

    write(*,*) ''
    write(*,*) 'Printing input values...'

    ! ~~~~~~~~~~~ simulation params ~~~~~~~~~~~~~~~~~~~
    write(*,*) ''
    write(*,*) 'Simulation parameters:'
    write(*,'(a32, i4, a10)') str_padding('OMP running with',31), &
                                    omp_get_max_threads(), str_padding('threads',9)
    write(*,'(a32, f9.4, a3)') str_padding('Estimated memory usage',31), &
                                   24.0*28.0*nPart/1024.0**2.0, 'MB'
    write(*,*) ''
    write(*,'(a24, i20)') str_padding('PRNG initializing seed',23), seed
    write(*,'(a32, i12)') str_padding('Number of atoms',31), nPart
    write(*,'(a26, f18.8)') str_padding('Initial box length (angs)',25), L
    write(*,'(a26, f18.8)') str_padding('Temperature (kelvin)',25), temp
    write(*,'(a26, e18.8)') str_padding('Pressure (bar)',25), pres
    write(*,'(a32, i12)') str_padding('MCMC walkers',31), walkers
    write(*,'(a32, i12)') str_padding('MCMC burn in steps',31), burn
    write(*,'(a32, i12)') str_padding('MCMC production steps',31), prod
    write(*,'(a22, i12, a10)') str_padding('Saving ',21), snaps, str_padding('snapshots',9)
    write(*,*) ''

    end subroutine init

! set up variables
subroutine setup_vars()
    implicit none
    allocate(x(nPart,3))
    write(*,*) 'Initializing PRNG...'
    call init_PRNG(seed)
    write(*,*) 'Success!'
    write(*,*) ''
    pres = pres/evang_to_bar ! convert bar to meV/angs^3

    if (trim(inimode)=='resume') then
        write(*,'(a50)') str_padding('Reading initial positions from ' // fstate, 49)
        write(*,*) ''
        open(unit=332, file=fstate)
        do n=1,8
            read(332,*)
        end do
        do i=1,nPart
            read(332,*) dummy, x(i,:)
        end do
        close(332)
    else if (trim(inimode)=='read') then
        write(*,'(a50)') str_padding('Reading initial positions from ' // fpos, 49)
        write(*,*) ''
        open(unit=332, file=fpos)
        read(332,*)
        read(332,*)
        do i=1,nPart
            read(332,*) dummy, x(i,:)
        end do
        close(332)
    else if (trim(inimode)=='random') then
        write(*,'(a50)') str_padding('Generating initial positions at random', 49)
        write(*,*) ''
        call random_number(x)
        x = x*L
    end if

    call ini_pbc_table()
    energy = compute_poten(nPart, L, x)

    end subroutine setup_vars

! write burn in step output
subroutine bin_step()
    write(*,'(a5, i6)') str_padding('Step',4), n
    write(*,*) ''
    write(*,'(a32, e16.8)') str_padding('Energy',31), energy/nPart
    write(*,'(a32, f16.8)') str_padding('BoxLength',31), L
    write(*,'(a32, f16.8)') str_padding('Density',31), density
    write(*,'(a32, e16.8)') str_padding('VirialPress',31), vpress
    write(*,*) ''
    write(*,*) '----------------------'
    write(*,*) ''
    end subroutine bin_step

! save final step
subroutine save_fstep()
    open(file='finalstate.out', unit=237)
    write(237,*) ''
    write(237,'(a32, e16.8)') str_padding('Energy',31), energy/nPart
    write(237,'(a32, f16.8)') str_padding('BoxLength',31), L
    write(237,'(a32, f16.8)') str_padding('Density',31), density
    write(237,'(a32, e16.8)') str_padding('VirialPress',31), vpress
    write(237,*) ''
    write(237,*) 'Positions'
    write(237,*) ''
    do i=1,nPart
        write(237,'(a4, 3f16.8)') asp, x(i,:)
    end do
    close(237)
    end subroutine save_fstep

! precompute pbc vectors table
subroutine ini_pbc_table()
    implicit none
    integer :: kx, ky, kz
    i = 0
    do kx=-1,1
    do ky=-1,1
    do kz=-1,1
    if (kx.ne.0 .or. ky.ne.0 .or. kz.ne.0) then
        i = i + 1
        pbc_table(i,:) = (/kx, ky, kz/)
    end if
    end do
    end do
    end do
    end subroutine ini_pbc_table
    
end program main
