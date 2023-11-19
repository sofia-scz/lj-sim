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
    integer :: i, j, p, ptot, npm, npp
    real(dp) :: q, d2, s2, rc2, rij(3), vpress
    vpress = zero
    rc2 = (sigma*rcut)**2.0_xp
    s2 = sigma**2.0_xp

    ptot = npart*(npart-1)/2-1
    npm = npart - 1
    npp = npart + 1

    !$omp parallel do private(p,i,j,d2,q,rij) reduction(+:vpress)
    do p=0,ptot
        ! get particle indices
        i = p/npm + 1
        j = mod(p,npm) + 1
        if (j.lt.i) then
            i = npp - i
            j = npp - j
        else 
            j = j + 1
        end if
        ! perform computations
        rij = x(i,:) - x(j,:)
        rij = rij - L*int(2.0_dp*rij/L)
        d2 = dot_product(rij, rij)
        if (d2 .lt. rc2) then
            q = (s2/d2)*(s2/d2)*(s2/d2)
            vpress = vpress + (2.0_dp*q*q - q)
        end if
    end do

    vpress = (8.0_dp*e0*vpress + npart*temp*kB)/L**3.0_xp

    end function compute_vpress


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function compute_poten(L, x) result(erg)
        real(dp), intent(in) :: L, x(npart,3)
        real(dp) :: erg
        integer :: i, j, p, ptot, npm, npp
        real(dp) :: q, d2, s2, rc2, rij(3), ercut
        erg = zero
        rc2 = (sigma*rcut)**2.0_xp
        s2 = sigma**2.0_xp
        ercut = (s2/rc2)**6.0_xp - (s2/rc2)**3.0_xp

        ptot = npart*(npart-1)/2-1
        npm = npart - 1
        npp = npart + 1
    
        !$omp parallel do private(p,i,j,d2,q,rij) reduction(+:erg)
        do p=0,ptot
            ! get particle indices
            i = p/npm + 1
            j = mod(p,npm) + 1
            if (j.lt.i) then
                i = npp - i
                j = npp - j
            else 
                j = j + 1
            end if
            ! perform computations
            rij = x(i,:) - x(j,:)
            rij = rij - L*int(2.0_dp*rij/L)
            d2 = dot_product(rij, rij)
            if (d2 .lt. rc2) then
                q = (s2/d2)*(s2/d2)*(s2/d2)
                erg = erg + (q*q - q - ercut)
            end if
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
            q = (s2/d2)*(s2/d2)*(s2/d2)
            de = de + (q - q*q + ercut)
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
            q = (s2/d2)*(s2/d2)*(s2/d2)
            de = de + (q*q - q - ercut)
        end if
    end if
    end do

    ! finish
    de = 4.0_dp*e0*de
    end function erg_diff

end module potential
