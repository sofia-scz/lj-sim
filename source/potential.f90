module potential
use global, only: dp, xp, zero, kB, &
                  npart, sigma, e0, rcut
use omp_lib
implicit none
private
public compute_poten, erg_diff, compute_vpress
contains

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function compute_vpress(temp, L, x) result(vpress)
    real(dp), intent(in) :: temp, L, x(npart,3)
    integer :: i, j
    real(dp) :: q, d2, s2, rc2, rij(3), vpress
    vpress = zero
    rc2 = (sigma*rcut)**2.0_xp
    s2 = sigma**2.0_xp

    !$omp parallel do private(i,j,d2,q,rij) reduction(+:vpress)
    do i=1,npart
    do j=i+1,npart
        rij = x(i,:) - x(j,:)
        rij = rij - L*int(2.0_dp*rij/L)
        d2 = dot_product(rij, rij)
        if (d2 .lt. rc2) then
            q = (s2/d2)**3.0_xp
            vpress = vpress + (2.0_dp*q**2.0_xp - q)
        end if
    end do
    end do

    vpress = (8.0_dp*e0*vpress + npart*temp*kB)/L**3.0_xp

    end function compute_vpress


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function compute_poten(L, x) result(erg)
        real(dp), intent(in) :: L, x(npart,3)
        real(dp) :: erg
        integer :: i, j
        real(dp) :: q, d2, s2, rc2, rij(3), ercut
        erg = zero
        rc2 = (sigma*rcut)**2.0_xp
        s2 = sigma**2.0_xp
        ercut = (s2/rc2)**6.0_xp - (s2/rc2)**3.0_xp
    
        !$omp parallel do private(i,j,d2,q,rij) reduction(+:erg)
        do i=1,npart
        do j=i+1,npart
            rij = x(i,:) - x(j,:)
            rij = rij - L*int(2.0_dp*rij/L)
            d2 = dot_product(rij, rij)
            if (d2 .lt. rc2) then
                q = (s2/d2)**3.0_xp
                erg = erg + (q**2.0_xp - q - ercut)
            end if
        end do
        end do

        erg = 4.0_dp*e0*erg
    
        end function compute_poten


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! erg perturbation
function erg_diff(i, L, x, xi_new) result(de)
    integer, intent(in) :: i
    real(dp), intent(in) :: L, x(npart,3), xi_new(3)
    real(dp) :: de
    integer :: j
    real(dp) :: q, d2, s2, rc2, rij(3), ercut
    de = zero
    rc2 = (sigma*rcut)**2.0_xp
    s2 = sigma**2.0_xp
    ercut = (s2/rc2)**6.0_xp - (s2/rc2)**3.0_xp    

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
            q = (s2/d2)**3.0_xp
            de = de + (q - q**2.0_xp + ercut)
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
            q = (s2/d2)**3.0_xp
            de = de + (q**2.0_xp - q - ercut)
        end if
    end if
    end do

    ! finish
    de = 4.0_dp*e0*de
    end function erg_diff

end module potential
