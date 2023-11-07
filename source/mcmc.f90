module mcmc
use global, only: dp, xp, zero
use random, only: randint
use potential, only: compute_poten, erg_diff
implicit none
real(dp), parameter :: xvar=0.02_dp, vvar=0.001_dp
private
public mcmc_move_NPT
contains

subroutine mcmc_move_NPT(nPart, temp, pres, L, x, energy)
    implicit none
    integer, intent(in) :: nPart
    real(dp), intent(in) :: temp, pres
    real(dp), intent(inout) :: L, x(nPart,3), energy
    integer :: i, j, k
    real(dp) :: u, e0, en, xj_old(3), xj_new(3), dx(3), &
                 v0, vn, Ln, logv

    e0 = energy

    ! particle positions perturbations
    do i=1,nPart
        ! select random particle
        j = randint(nPart) + 1

        ! save its old position
        xj_old = x(j,:)

        ! generate new position
        call random_number(dx)
        dx = (dx - 0.5_dp)*xvar*L
        xj_new = x(j,:) + dx
        
        ! fix particles escaping the box
        do k=1,3
            if (xj_new(k) .lt. zero) then
                xj_new(k) = xj_new(k) + L
            else if (xj_new(k) .ge. L) then
                xj_new(k) = xj_new(k) - L
            end if
        end do
        
        ! compute new erg
        en = e0 + erg_diff(nPart, j, L, x, xj_new)

        ! acceptance test
        call random_number(u)
        if (u .lt. exp((e0-en)/temp)) then
            e0 = en
            x(j,:) = xj_new
        end if
    end do

    energy = e0

    ! volume perturbations
    v0 = L**3.0_xp
    call random_number(u)

    logv = log(v0) + (u-0.5_dp)*vvar
    vn = exp(logv)
    Ln = vn**(1.0_xp/3.0_xp)

    en = compute_poten(nPart, Ln, x/L*Ln)
    call random_number(u)

    if (u .lt. exp(-(en - e0 + pres*(vn-v0))/temp &
                    -(nPart+1)*log(vn/v0))) then
        x = x/L*Ln
        L = Ln
        energy = en
    end if

    end subroutine mcmc_move_NPT

end module mcmc