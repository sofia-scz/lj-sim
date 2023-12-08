program ljmd
use global, only:       dp, int64, & 
                        kB, amuangs_to_kglt, evangs_to_bar, &
                        amu_speed, &
                        npart, mass, e0, sigma, rcut
use random, only:       init_PRNG, randn_matrix
use io, only:           md_readin, readsys
use dynamics, only:     NVE_verlet, langevin
use potential, only:    get_poten, get_force, get_press
use omp_lib
implicit none

real(dp), allocatable :: x(:,:), v(:,:), a(:,:), f(:,:)
real(dp) :: L, dt, temp, press, target_temp, target_press, gamma, &
            toterg, poterg, kinerg, t0, tf
integer :: i, n, relax, steps, snaps
integer(int64) :: seed
character(len=2) :: asp, dummy
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

do n=1,relax
    v = 0
    call langevin(dt, L, target_temp, gamma, x, v, a, f)
end do

! open data dump files
open(file='table.out', unit=236)
write(236,'(a81)') 'Energy          Potential         Kinetic        &
                    &Temperature       Pressure'
open(unit=145, file='snapshots.out')

do n=1,steps
    call langevin(dt, L, target_temp, gamma, x, v, a, f)
    call write_table()
    ! write snapshots
    if (mod(n-mod(steps,snaps),snaps).eq.0) then
        call write_snap()
        call flush()
    end if
    ! force data dump every 1%
    if (mod(n-mod(steps,int(steps/100)),int(steps/100)).eq.0) then
        call flush()
    end if
end do

close(145)
close(236)
tf = omp_get_wtime()
write(*,'(a33,x,f12.4,a4)') 'Simulation finished! Time spent:', tf-t0, 'sec'
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
    call md_readin(seed, L, dt, target_temp, target_press, gamma, relax, steps, snaps, inimode)
    call readsys(npart, mass, e0, sigma, rcut, asp)
    write(*,*) 'Done.'

    write(*,*)
    write(*,*) 'Printing input values...'

    ! ~~~~~~~~~~~ simulation params ~~~~~~~~~~~~~~~~~~~
    write(*,*)
    write(*,*) 'Simulation parameters:'
    write(*,'(a20,18x,i4,a8)') 'OpenMP running with', omp_get_max_threads(), 'threads'
    write(*,'(a23,16x,f8.3,a3)') 'Estimated memory usage', 3*npart*4/1024.0**2*1.2, 'MB'
    write(*,*)
    write(*,'(a15,31x,a4)') 'Atomic species', asp
    write(*,'(a18,20x,f12.4)') 'Atomic mass (amu)', mass
    write(*,'(a16,22x,i12)') 'Number of atoms', npart
    write(*,'(a22,16x,g12.6)') 'Lennard Jones e0 (eV)', e0
    write(*,'(a27,11x,f12.4)') 'Lennard Jones sigma (angs)', sigma
    write(*,'(a26,16x,f8.4)') 'Cutoff radius (LJ sigmas)', rcut
    write(*,*)
    write(*,'(a23,7x,i20)') 'PRNG initializing seed', seed
    write(*,'(a26,12x,f12.4)') 'Initial box length (angs)', L
    write(*,'(a28,10x,f12.4)') 'Target temperature (kelvin)', target_temp
    write(*,'(a22,16x,f12.4)') 'Target pressure (bar)', target_press
    write(*,'(a16,22x,f12.4)') 'Density (Kg/lt)', mass*amuangs_to_kglt*npart/L**3
    write(*,'(a15,23x,f12.4)') 'Langevin gamma', gamma
    write(*,*)
    write(*,'(a15,23x,f12.2)') 'Time step (fs)', dt
    write(*,'(a9,29x,i12)') 'MD steps', steps
    write(*,*) 'Saving snapshots every', snaps, 'steps'
    write(*,*)
    write(*,*)

    end subroutine init

! set up variables
subroutine setup_vars()
    implicit none
    ! allocate positions matrix
    allocate(x(npart,3))
    allocate(v(npart,3))
    allocate(a(npart,3))
    allocate(f(npart,3))
    mass = mass*amu_speed

    ! initialize PRNG from seed
    write(*,*) 'Initializing PRNG...'
    call init_PRNG(seed)
    write(*,*) 'Success!'
    write(*,*)

    ! get position values
    if (trim(inimode)=='read') then
        write(*,*) 'Reading initial conditions from ' // fpos
        write(*,*)
        open(unit=332, file=fpos)
        read(332,*)
        read(332,*)
        do i=1,npart
            read(332,*) dummy, x(i,:), v(i,:), f(i,:)
        end do
        close(332)
    else if (trim(inimode)=='random') then
        write(*,*) 'Generating initial conditions at random'
        write(*,*)
        call random_number(x)
        x = x*L
        v = randn_matrix(npart,3)*(3.0_dp*kB*target_temp/mass)**0.5_dp
        f = get_force(L, x)
    end if

    ! set up some variables
    a = f/mass
    poterg = get_poten(L, x)
    kinerg = get_kinetic()
    temp = get_temp(kinerg)
    end subroutine setup_vars

! write to table.out
subroutine write_table()
    implicit none

    ! compute
    poterg = get_poten(L, x)
    kinerg = get_kinetic()
    toterg = poterg + kinerg
    temp = get_temp(kinerg)
    press = get_press(temp, L, x)*evangs_to_bar

    ! write line
    write(236,'(e20.12,x,e20.12,x,e20.12,x,e20.12,x,e20.12)') &
        toterg, poterg, kinerg, temp, press
    end subroutine write_table

! write snapshot
subroutine write_snap()
    implicit none

    ! header
    write(145,*) npart
    write(145,'(a54,a10,f7.3,a13,f7.3,a13,f7.3,a13)') &
            'Properties=species:S:1:pos:R:3:momenta:R:3:forces:R:3 ', &
            'Lattice="', L, ' 0.0 0.0 0.0 ', L, ' 0.0 0.0 0.0 ', L, &
            '" pbc="T T T"'

    ! particles
    do i=1,npart
        write(145,'(a3,9f16.8)') asp, x(i,:), v(i,:), f(i,:)
    end do
    end subroutine write_snap

! compute kinetic energy
function get_kinetic() result(kin)
    implicit none
    real(dp) :: kin
    
    kin = 0.0d0
    do i=1,npart
        kin = kin + dot_product(v(i,:),v(i,:))
    end do
    kin = 0.5_dp*mass*kin
    end function get_kinetic

! compute temperature
function get_temp(kin) result(tcalc)
    implicit none
    real(dp), intent(in) :: kin
    real(dp) :: tcalc

    tcalc = 2.0_dp*kin/(3.0_dp*npart*kB)

    end function get_temp

end program ljmd
