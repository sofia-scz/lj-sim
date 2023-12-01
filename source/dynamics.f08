module dynamics
use global, only: dp, npart, mass, kB
use potential, only: get_force
use random, only: randn_matrix
private
public NVE_verlet, langevin
contains

subroutine NVE_verlet(dt, L, x, v, a, f)
    implicit none
    real(dp), intent(in) :: dt, L
    real(dp), dimension(npart,3), intent(inout) :: x, v, a, f

    v = v + (0.5_dp*dt)*a
    x = x + dt*v
    x = modulo(x,L)
    f = get_force(L, x)
    a = f/mass
    v = v + (0.5_dp*dt)*a

    end subroutine NVE_verlet

subroutine langevin(dt, L, temp, y, x, v, a, f)
    implicit none
    real(dp), intent(in) :: dt, L, temp, y
    real(dp), dimension(npart,3), intent(inout) :: x, v, a, f

    v = v + (0.5_dp*dt)*a
    x = x + dt*v
    x = modulo(x,L)
    f = get_force(L, x)
    a = f/mass - y*v + (2.0_dp*y*kB*temp/(mass*dt))**0.5_dp*randn_matrix(npart,3)
    v = v + (0.5_dp*dt)*a

    end subroutine langevin

end module dynamics