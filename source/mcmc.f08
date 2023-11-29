module mcmc
use global, only: dp, kB, npart, xvar, vvar
use random, only: randint
use potential, only: get_poten, erg_diff
implicit none
private
public NpT_step
contains

subroutine NpT_step(temp, target_press, L, x, energy, xtry, vtry, xacc, vacc)
    implicit none
    real(dp), intent(in) :: temp, target_press
    real(dp), intent(inout) :: L, x(npart,3), energy, &
                                xtry, vtry, xacc, vacc
    integer :: i, ri, accept

    do i=1,npart+1
        ri = randint(npart+1)
        if (ri.eq.0) then
            vtry = vtry + 1
            call volstretch_move(temp, target_press, L, x, energy, accept)
            vacc = vacc + accept
        else
            xtry = xtry + 1
            call displace_move(temp, L, x, energy, accept)
            xacc = xacc + accept
        end if
    end do

    end subroutine NpT_step

! attempt to displace 1 particle
subroutine displace_move(temp, L, x, energy, accept)
    implicit none
    real(dp), intent(in) :: temp, L
    real(dp), intent(inout) :: x(npart,3), energy
    integer, intent(out) :: accept
    integer :: i, k
    real(dp) :: u, eold, enew, xi_old(3), xi_new(3), dx(3), beta

    eold = energy
    beta = 1.0_dp/(kB*temp)
    accept = 0

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
        if (xi_new(k) .lt. 0.0_dp) then
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
        accept = 1
    end if

    end subroutine displace_move

! attempt to stretch the whole simulation box
subroutine volstretch_move(temp, target_press, L, x, energy, accept)
    implicit none
    real(dp), intent(in) :: temp, target_press
    real(dp), intent(inout) :: L, x(npart,3), energy
    integer, intent(out) :: accept
    real(dp) :: u, eold, enew, v0, vn, Ln, logv, beta

    eold = energy
    beta = 1.0_dp/(kB*temp)
    accept = 0

    ! volume perturbations
    v0 = L**3
    call random_number(u)

    logv = log(v0) + (u-0.5_dp)*2.0_dp*vvar
    vn = exp(logv)
    Ln = vn**(1.0_dp/3.0_dp)

    enew = get_poten(Ln, x/L*Ln)
    call random_number(u)

    if (u .lt. exp(-beta*(enew - eold + target_press*(vn-v0)) &
                        + npart*log(vn/v0))) then
        x = x/L*Ln
        L = Ln
        energy = enew
        accept = 1
    end if

    end subroutine volstretch_move

end module mcmc