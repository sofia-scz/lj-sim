module potential
use global, only: dp, kB, &
                  npart, sigma, e0, rcut
use omp_lib
implicit none
private
public get_poten, erg_diff, get_force, get_press
contains

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function get_press(temp, L, x) result(press)
    real(dp), intent(in) :: temp, L, x(npart,3)
    integer :: i, j
    real(dp) :: q, d2, s2, rc2, rij(3), press
    press = 0.0_dp
    rc2 = (sigma*rcut)**2
    s2 = sigma**2

    !$omp parallel do private(i,j,d2,q,rij) schedule(static,1) reduction(+:press)
    do i=1,npart-1
    do j=i+1,npart
        ! perform computations
        rij = x(i,:) - x(j,:)
        rij = rij - L*int(2.0_dp*rij/L)
        d2 = dot_product(rij, rij)
        if (d2 .lt. rc2) then
            q = (s2/d2)**3
            press = press + (2.0_dp*q**2 - q)
        end if
    end do
    end do

    press = (8.0_dp*e0*press + npart*temp*kB)/L**3

    end function get_press


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function get_force(L, x) result(force)
    real(dp), intent(in) :: L, x(npart,3)
    integer :: i, j
    real(dp) :: q, d2, s2, rc2, rij(3), force(npart,3), fij(3)
    force = 0.0_dp
    rc2 = (sigma*rcut)**2
    s2 = sigma**2

    !$omp parallel do private(i,j,d2,q,rij,fij) schedule(static,1) reduction(+:force)
    do i=1,npart-1
    do j=i+1,npart
        ! perform computations
        rij = x(i,:) - x(j,:)
        rij = rij - L*int(2.0_dp*rij/L)
        d2 = dot_product(rij, rij)
        if (d2 .lt. rc2) then
            q = (s2/d2)**3
            fij = rij/d2*(2.0_dp*q**2 - q)
            force(i,:) = force(i,:) + fij
            force(j,:) = force(j,:) + (-fij)
        end if
    end do
    end do

    force = 24.0_dp*e0*force

    end function get_force


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function get_poten(L, x) result(erg)
        real(dp), intent(in) :: L, x(npart,3)
        real(dp) :: erg
        integer :: i, j
        real(dp) :: q, d2, s2, rc2, rij(3), ercut
        erg = 0.0_dp
        rc2 = (sigma*rcut)**2
        s2 = sigma**2
        ercut = (s2/rc2)**6 - (s2/rc2)**3
    
        !$omp parallel do private(i,j,d2,q,rij) schedule(static,1) reduction(+:erg)
        do i=1,npart-1
        do j=i+1,npart
            ! perform computations
            rij = x(i,:) - x(j,:)
            rij = rij - L*int(2.0_dp*rij/L)
            d2 = dot_product(rij, rij)
            if (d2 .lt. rc2) then
                q = (s2/d2)**3
                erg = erg + (q**2 - q - ercut)
            end if
        end do
        end do

        erg = 4.0_dp*e0*erg
    
        end function get_poten


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! erg perturbation
function erg_diff(i, L, x, xi_new) result(de)
    integer, intent(in) :: i
    real(dp), intent(in) :: L, x(npart,3), xi_new(3)
    real(dp) :: de
    integer :: j
    real(dp) :: q, d2, s2, rc2, rij(3), ercut
    de = 0.0_dp
    rc2 = (sigma*rcut)**2
    s2 = sigma**2
    ercut = (s2/rc2)**6 - (s2/rc2)**3    

    !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! substract old interactions
    !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !$omp parallel do private(j,d2,q,rij) reduction(+:de)
    do j=1,npart
    if (j.ne.i) then
        rij = x(i,:) - x(j,:)
        rij = rij - L*int(2.0_dp*rij/L)
        d2 = dot_product(rij, rij)
        if (d2 .lt. rc2) then
            q = (s2/d2)**3
            de = de + (q - q**2 + ercut)
        end if
    end if
    end do

    !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! add new interactions
    !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !$omp parallel do private(j,d2,q,rij) reduction(+:de)
    do j=1,npart
    if (j.ne.i) then
        rij = xi_new - x(j,:)
        rij = rij - L*int(2.0_dp*rij/L)
        d2 = dot_product(rij, rij)
        if (d2 .lt. rc2) then
            q = (s2/d2)**3
            de = de + (q**2 - q - ercut)
        end if
    end if
    end do

    ! finish
    de = 4.0_dp*e0*de
    end function erg_diff

end module potential
