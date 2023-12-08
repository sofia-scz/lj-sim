program ljmc
use global, only:       dp, int64, & 
                        amuangs_to_kglt, evangs_to_bar, &
                        npart, mass, e0, sigma, rcut, &
                        xvar, vvar
use random, only:       init_PRNG
use io, only:           mc_readin, readsys
use mcmc, only:         NpT_step
use potential, only:    get_poten, get_press
use omp_lib
implicit none

real(dp), allocatable :: x(:,:)
real(dp) :: L, target_temp, target_press, t0, tf, energy, & 
            density, press, xtry, vtry, xacc, vacc
integer :: i, n, burn, prod, snaps
integer(int64) :: seed
character(len=2) :: asp, dummy
character(len=9) :: number
character(len=10) :: inimode
character(len=14) :: fpos='positions.in'

! call init
call init()
call flush()

! set up variables
call setup_vars()
call flush()

! start computations
write(*,*)
write(*,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(*,*) '~~~~~~~~~~~~ Starting computations ~~~~~~~~~~~~'
write(*,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(*,*)
call flush()

t0 = omp_get_wtime()
write(*,*) 'Starting burn in...'
call flush()

do n=1,burn
    call NpT_step(target_temp, target_press, L, x, energy, xtry, vtry, xacc, vacc)
end do

tf = omp_get_wtime()
write(*,'(a30,4x,f12.4,a4)') 'Burn in finished! Time spent:', tf-t0, 'sec'
call flush()

t0 = omp_get_wtime()
write(*,*)
write(*,*) 'Starting production...'
call flush()

! open data dump files
open(file='table.out', unit=236)
write(236,'(a82)') 'Energy        BoxLength        Density        &
                    &Pressure        XACR        VACR'
open(unit=145, file='snapshots.out')

do n=1,prod
    call NpT_step(target_temp, target_press, L, x, energy, xtry, vtry, xacc, vacc)
    call write_table()
    ! write snapshots
    if (mod(n-mod(prod,int(prod/snaps)),int(prod/snaps)).eq.0) then
        call write_snap()
        call flush()
    end if
    ! force data dump every 1%
    if (mod(n-mod(prod,int(prod/100)),int(prod/100)).eq.0) then
        call flush()
    end if
end do

close(145)
close(236)
tf = omp_get_wtime()
write(*,'(a33,x,f12.4,a4)') 'Production finished! Time spent:', tf-t0, 'sec'
write(*,*)

write(*,*)
write(*,*) '      ~(˘▾˘~) Execution completed (~˘▾˘)~       '
write(*,*)
write(*,*) 'Have a nice day!'
write(*,*)

contains

! initialize program
subroutine init()
    implicit none
    ! start writing output
    write(*,*)
    write(*,*) '-------------------------------------------------------------'
    write(*,*) '| Welcome to my Lennard Jones molecular simulation!!        |'
    write(*,*) '|                                                           |'
    write(*,*) '| This program was written by C. Dacal and S. Scozziero for |'
    write(*,*) '| the course "Introduccion a la Simulacion Computacional"   |'
    write(*,*) '| taught by Dr. Pastorino at Instituto Sabato.              |'
    write(*,*) '|                                                           |'
    write(*,*) '| We hope you enjoy using our program (✿◠‿◠)                |'
    write(*,*) '-------------------------------------------------------------'

    write(*,*)
    write(*,*) 'Reading input files...'
    call mc_readin(seed, L, target_temp, target_press, xvar, vvar, burn, prod, snaps, inimode)
    call readsys(npart, mass, e0, sigma, rcut, asp)
    write(*,*) 'Done.'

    write(*,*)
    write(*,*) 'Printing input values...'

    ! ~~~~~~~~~~~ simulation params ~~~~~~~~~~~~~~~~~~~
    write(*,*)
    write(*,*) 'Simulation parameters:'
    write(*,'(a20,18x,i4,a8)') 'OpenMP running with', omp_get_max_threads(), 'threads'
    write(*,'(a23,16x,f8.3,a3)') 'Estimated memory usage', 3*npart/1024.0**2*1.2, 'MB'
    write(*,*)
    write(*,'(a15,31x,a4)') 'Atomic species', asp
    write(*,'(a18,20x,f12.4)') 'Atomic mass (amu)', mass
    write(*,'(a16,22x,i12)') 'Number of atoms', npart
    write(*,'(a22,16x,g12.6)') 'Lennard Jones e0 (eV)', e0
    write(*,'(a27,11x,f12.6)') 'Lennard Jones sigma (angs)', sigma
    write(*,'(a26,16x,f8.4)') 'Cutoff radius (LJ sigmas)', rcut
    write(*,*)
    write(*,'(a23,7x,i20)') 'PRNG initializing seed', seed
    write(*,'(a26,12x,f12.6)') 'Initial box length (angs)', L
    write(*,'(a28,10x,f12.6)') 'Target Temperature (kelvin)', target_temp
    write(*,'(a22,16x,f12.6)') 'Target Pressure (bar)', target_press
    write(*,*)
    write(number,'(g9.3)') xvar*100
    write(*,*) 'Max Particle Displacement (% of box length) ' // number
    write(number,'(g9.3)') vvar*100
    write(*,*) 'Max Volume Stretch (% of box volume)        ' // number
    write(*,'(a19,19x,i12)') 'MCMC burn in steps', burn
    write(*,'(a22,16x,i12)') 'MCMC production steps', prod
    write(*,*) 'Saving snapshots every', snaps, 'steps'
    write(*,*)

    end subroutine init

! set up variables
subroutine setup_vars()
    implicit none
    ! allocate positions matrix
    allocate(x(npart,3))

    ! initialize PRNG from seed
    write(*,*) 'Initializing PRNG...'
    call init_PRNG(seed)
    write(*,*) 'Success!'
    write(*,*)

    ! get position values
    if (trim(inimode)=='read') then
        write(*,*) 'Reading initial positions from ' // fpos
        write(*,*)
        open(unit=332, file=fpos)
        read(332,*)
        read(332,*)
        do i=1,npart
            read(332,*) dummy, x(i,:)
        end do
        close(332)
    else if (trim(inimode)=='random') then
        write(*,*) 'Generating initial positions at random'
        write(*,*)
        call random_number(x)
        x = x*L
    end if

    ! set up some variables
    energy = get_poten(L, x)
    target_press = target_press/evangs_to_bar
    xacc = 0
    vacc = 0
    xtry = 1
    vtry = 1

    end subroutine setup_vars

! write to table.out
subroutine write_table()
    implicit none

    ! compute
    density = npart*mass/L**3*amuangs_to_kglt
    press = get_press(target_temp, L, x)*evangs_to_bar

    ! write line
    write(236,'(e20.12,x,e20.12,x,e20.12,x,e20.12,x,f10.6,x,f10.6)') &
        energy/npart*1000, L, density, press, xacc/xtry, vacc/vtry
    end subroutine write_table

! write write_snapshot
subroutine write_snap()
    implicit none

    ! header
    write(145,*) npart
    write(145,'(a10,f7.3,a13,f7.3,a13,f7.3,a13)') 'Lattice="', L, &
            ' 0.0 0.0 0.0 ', L, ' 0.0 0.0 0.0 ', L, '" pbc="T T T"'

    ! particles
    do i=1,npart
        write(145,'(a3,x,f16.8,x,f16.8,x,f16.8)') asp, x(i,:)
    end do
    end subroutine write_snap
    
end program ljmc
