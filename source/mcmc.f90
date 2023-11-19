module mcmc
use global, only: dp, xp, zero, kB, npart, xvar, vvar
use random, only: randint
use potential, only: compute_poten, erg_diff
implicit none
private
public NpT_step
contains

subroutine NpT_step(temp, pres, L, x, energy, xtry, vtry, xacc, vacc)
    implicit none
    real(dp), intent(in) :: temp, pres
    real(dp), intent(inout) :: L, x(npart,3), energy, &
                                xtry, vtry, xacc, vacc
    logical :: accept
    integer :: i, ri

    do i=1,npart+1
        ri = randint(npart+1)
        if (ri.eq.0) then
            vtry = vtry + 1
            call volstretch_move(temp, pres, L, x, energy, accept)
            if (accept) then
                vacc = vacc + 1
            end if
        else
            xtry = xtry + 1
            call displace_move(temp, L, x, energy, accept)
            if (accept) then
                xacc = xacc + 1
            end if
        end if
    end do

    end subroutine NpT_step

! attempt to displace 1 particle
subroutine displace_move(temp, L, x, energy, accept)
    implicit none
    real(dp), intent(in) :: temp, L
    real(dp), intent(inout) :: x(npart,3), energy
    logical, intent(out) :: accept
    integer :: i, k
    real(dp) :: u, eold, enew, xi_old(3), xi_new(3), dx(3), beta

    eold = energy
    beta = 1.0_dp/(kB*temp)
    accept = .false.

    ! select random particle
    i = randint(npart) + 1

    ! save its old position
    xi_old = x(i,:)

    ! generate new position
    call random_number(dx)
    dx = (dx - 0.5_dp)*xvar*L
    xi_new = x(i,:) + dx
    
    ! fix particles escaping the box
    do k=1,3
        if (xi_new(k) .lt. zero) then
            xi_new(k) = xi_new(k) + L
        else if (xi_new(k) .ge. L) then
            xi_new(k) = xi_new(k) - L
        end if
    end do
    
    ! compute new erg
    enew = eold + erg_diff(i, L, x, xi_new)

    ! acceptance test
    call random_number(u)
    if (u .lt. exp((eold-enew)*beta)) then
        energy = enew
        x(i,:) = xi_new
        accept = .true.
    end if

    end subroutine displace_move

! attempt to stretch the whole simulation box
subroutine volstretch_move(temp, pres, L, x, energy, accept)
    implicit none
    real(dp), intent(in) :: temp, pres
    real(dp), intent(inout) :: L, x(npart,3), energy
    logical, intent(out) :: accept
    real(dp) :: u, eold, enew, v0, vn, Ln, logv, beta

    eold = energy
    beta = 1.0_dp/(kB*temp)
    accept = .false.

    ! volume perturbations
    v0 = L**3.0_xp
    call random_number(u)

    logv = log(v0) + (u-0.5_dp)*2.0_dp*vvar
    vn = exp(logv)
    Ln = vn**(1.0_xp/3.0_xp)

    enew = compute_poten(Ln, x/L*Ln)
    call random_number(u)

    if (u .lt. exp(-beta*(enew - eold + pres*(vn-v0)) & 
                        + npart*log(vn/v0))) then
        x = x/L*Ln
        L = Ln
        energy = enew
        accept = .true.
    end if

    end subroutine volstretch_move

end module mcmc