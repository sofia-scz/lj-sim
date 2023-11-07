module potential
use global, only: dp, xp, int64, zero, kB
use omp_lib
implicit none
real(dp), parameter :: s_LJ=2.76125_dp, &      ! angs
                       e0_LJ=3.63242_dp, &     ! meV
                       mass=20.180             ! amu
integer :: pbc_table(26,3)
private
public compute_poten, erg_diff, compute_vpress, pbc_table, mass
contains

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function compute_force(nPart, L, x) result(force)
    integer, intent(in) :: nPart
    real(dp), intent(in) :: L, x(nPart,3)
    real(dp) :: force(nPart,3)

    integer :: i, j, k
    real(dp) :: q, d, rcut, x_PBC(nPart*26,3), f(3), rij(3)
    force = zero
    rcut = L*0.99_dp

    ! make periodic images
    do k=1,26
    do i=1,nPart
        x_PBC(i+(k-1)*nPart,:) = x(i,:) + L*pbc_table(k,:)
    end do
    end do 

    !$omp parallel do private(i,j,d,q,rij,f) reduction(+:force)
    do i=1,nPart
        do j=i+1,nPart
            rij = x(i,:) - x(j,:)
            d = norm2(rij)
            if (d .lt. rcut) then
                q = (s_LJ/d)**6.0_xp
                f = rij/d**2.0_xp*(2.0_dp*q**2.0_xp - q)
                force(i,:) = force(i,:) - f
                force(j,:) = force(j,:) + f
            end if
        end do
        do j=1,nPart*26
            rij = x(i,:) - x_PBC(j,:)
            d = norm2(rij)
            if (d .lt. rcut) then
                q = (s_LJ/d)**6.0_xp
                f = rij/d**2.0_xp*(2.0_dp*q**2.0_xp - q)
                force(i,:) = force(i,:) - f
            end if
        end do
    end do

    force = 24.0_dp*e0_LJ*force

    end function compute_force


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function compute_vpress(nPart, temp, L, x) result(vpress)
    integer, intent(in) :: nPart
    real(dp), intent(in) :: temp, L, x(nPart,3)
    integer :: i, j, k
    real(dp) :: q, d, rcut, x_PBC(nPart*26,3), rij(3), vpress
    vpress = zero
    rcut = L*0.99_dp

    ! make periodic images
    do k=1,26
    do i=1,nPart
        x_PBC(i+(k-1)*nPart,:) = x(i,:) + L*pbc_table(k,:)
    end do
    end do 

    !$omp parallel do private(i,j,d,q,rij) reduction(+:vpress)
    do i=1,nPart
        do j=i+1,nPart
            rij = x(i,:) - x(j,:)
            d = norm2(rij)
            if (d .lt. rcut) then
                q = (s_LJ/d)**6.0_xp
                vpress = vpress + (2.0_dp*q**2.0_xp - q)
            end if
        end do
        do j=1,nPart*26
            rij = x(i,:) - x_PBC(j,:)
            d = norm2(rij)
            if (d .lt. rcut) then
                q = (s_LJ/d)**6.0_xp
                vpress = vpress + (2.0_dp*q**2.0_xp - q)
            end if
        end do
    end do

    vpress = (8.0_dp*e0_LJ*vpress + nPart*temp*kB)/L**3.0_xp

    end function compute_vpress


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function compute_poten(nPart, L, x) result(erg)
        integer, intent(in) :: nPart
        real(dp), intent(in) :: L, x(nPart,3)
        real(dp) :: erg
    
        integer :: i, j, k
        real(dp) :: q, d, rcut, x_PBC(nPart*26,3)
        erg = zero
        rcut = L*0.99_dp
    
        ! make periodic images
        do k=1,26
        do i=1,nPart
            x_PBC(i+(k-1)*nPart,:) = x(i,:) + L*pbc_table(k,:)
        end do
        end do 
    
        !$omp parallel do private(i,j,d,q) reduction(+:erg)
        do i=1,nPart
            do j=i+1,nPart
                d = norm2(x(i,:) - x(j,:))
                if (d .lt. rcut) then
                    q = (s_LJ/d)**6.0_xp
                    erg = erg + (q**2.0_xp - q)
                end if
            end do
            do j=1,nPart*26
                d = norm2(x(i,:) - x_PBC(j,:))
                if (d .lt. rcut) then
                    q = (s_LJ/d)**6.0_xp
                    erg = erg + (q**2.0_xp - q)
                end if
            end do
        end do
    
        erg = 4.0_dp*e0_LJ*erg
    
        end function compute_poten


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! erg perturbation
function erg_diff(nPart, i, L, x, xi_new) result(de)
    integer, intent(in) :: nPart, i
    real(dp), intent(in) :: L, x(nPart,3), xi_new(3)
    real(dp) :: de

    integer :: j, k
    real(dp) :: q, d, rcut, x_PBC(nPart*26,3)
    de = zero
    rcut = L*0.99_dp

    !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! substract old interactions
    !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! make periodic images
    do k=1,26
    do j=1,nPart
        x_PBC(j+(k-1)*nPart,:) = x(j,:) + L*pbc_table(k,:)
    end do
    end do

    ! interactions between i and j particles
    do j=1,i-1
        d = norm2(x(i,:) - x(j,:))
        if (d .lt. rcut) then
            q = (s_LJ/d)**6.0_xp
            de = de - (q**2.0_xp - q)
        end if
    end do
    do j=i+1,nPart
        d = norm2(x(i,:) - x(j,:))
        if (d .lt. rcut) then
            q = (s_LJ/d)**6.0_xp
            de = de - (q**2.0_xp - q)
        end if
    end do

    ! interactions between i and j images
    !$omp parallel do private(j,d,q) reduction(+:de)
    do j=1,nPart*26
        d = norm2(x(i,:) - x_PBC(j,:))
        if (d .lt. rcut) then
            q = (s_LJ/d)**6.0_xp
            de = de + 2.0_dp*(-q**2.0_xp + q)
        end if
    end do
    
    !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! add new interactions
    !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! fix i periodic img
    do k=1,26
        x_PBC(i+(k-1)*nPart,:) = xi_new + L*pbc_table(k,:)
    end do 

    ! interactions between i and j particles
    do j=1,i-1
        d = norm2(xi_new - x(j,:))
        if (d .lt. rcut) then
            q = (s_LJ/d)**6.0_xp
            de = de + (q**2.0_xp - q)
        end if
    end do
    do j=i+1,nPart
        d = norm2(xi_new - x(j,:))
        if (d .lt. rcut) then
            q = (s_LJ/d)**6.0_xp
            de = de + (q**2.0_xp - q)
        end if
    end do

    ! interactions between i and j images
    !$omp parallel do private(j,d,q) reduction(+:de)
    do j=1,npart*26
        d = norm2(xi_new - x_PBC(j,:))
        if (d .lt. rcut) then
            q = (s_LJ/d)**6.0_xp
            de = de + 2.0_dp*(q**2.0_xp - q)
        end if
    end do

    ! finish
    de = 4.0_dp*e0_LJ*de
    end function erg_diff

! 
function test_cutoff(rij, rc) result(d)
    implicit none
    real(dp), intent(in) :: rij(3), rc
    real(dp) :: d

    if (abs(rij(1)).ge.rc .or. abs(rij(2)).ge.rc .or. abs(rij(3)).ge.rc) then
        d = -1.0
        return
    else
        d = norm2(rij)
        if (d .ge. rc) then
            d = -1.0
        end if
        return
    end if

    end function test_cutoff

end module potential

! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! !                     OLD
! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ! erg perturbation
! function erg_diff(nPart, i, L, x, xi_new) result(de)
!     integer, intent(in) :: nPart, i
!     real(dp), intent(in) :: L, x(nPart,3), xi_new(3)
!     real(dp) :: de

!     integer :: ii, j, k
!     real(dp) :: q, d, rcut, pbc(3), x_PBC(nPart*26,3)
!     de = zero
!     rcut = L*0.99_dp

!     !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ! substract old interactions
!     !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     ! make periodic images
!     do k=1,26
!     do j=1,nPart
!         x_PBC(j+(k-1)*nPart,:) = x(j,:) + L*pbc_table(k,:)
!     end do
!     end do

!     ! interactions between i and j particles
!     do j=1,i-1
!         d = norm2(x(i,:) - x(j,:))
!         if (d .lt. rcut) then
!             q = (s_LJ/d)**6.0_xp
!             de = de - (q**2.0_xp - q)
!         end if
!     end do
!     do j=i+1,nPart
!         d = norm2(x(i,:) - x(j,:))
!         if (d .lt. rcut) then
!             q = (s_LJ/d)**6.0_xp
!             de = de - (q**2.0_xp - q)
!         end if
!     end do

!     ! interactions between i images and j particles
!     do ii=i,26*nPart,nPart
!         do j=1,i-1
!             d = norm2(x_PBC(ii,:) - x(j,:))
!             if (d .lt. rcut) then
!                 q = (s_LJ/d)**6.0_xp
!                 de = de - (q**2.0_xp - q)
!             end if
!         end do
!         do j=i+1,nPart
!             d = norm2(x_PBC(ii,:) - x(j,:))
!             if (d .lt. rcut) then
!                 q = (s_LJ/d)**6.0_xp
!                 de = de - (q**2.0_xp - q)
!             end if
!         end do
!     end do

!     ! interactions between i and j images
!     !$omp parallel do private(j,d,q) reduction(+:de)
!     do j=1,nPart*26
!         d = norm2(x(i,:) - x_PBC(j,:))
!         if (d .lt. rcut) then
!             q = (s_LJ/d)**6.0_xp
!             de = de + (-q**2.0_xp + q)
!         end if
!     end do
    
!     !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ! add new interactions
!     !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     ! fix i periodic img
!     do k=1,26
!         x_PBC(i+(k-1)*nPart,:) = xi_new + L*pbc_table(k,:)
!     end do 

!     ! interactions between i and j particles
!     do j=1,i-1
!         d = norm2(xi_new - x(j,:))
!         if (d .lt. rcut) then
!             q = (s_LJ/d)**6.0_xp
!             de = de + (q**2.0_xp - q)
!         end if
!     end do
!     do j=i+1,nPart
!         d = norm2(xi_new - x(j,:))
!         if (d .lt. rcut) then
!             q = (s_LJ/d)**6.0_xp
!             de = de + (q**2.0_xp - q)
!         end if
!     end do

!     ! interactions between i images and j particles
!     do ii=i,26*nPart,nPart
!         do j=1,i-1
!             d = norm2(x_PBC(ii,:) - x(j,:))
!             if (d .lt. rcut) then
!                 q = (s_LJ/d)**6.0_xp
!                 de = de + (q**2.0_xp - q)
!             end if
!         end do
!         do j=i+1,nPart
!             d = norm2(x_PBC(ii,:) - x(j,:))
!             if (d .lt. rcut) then
!                 q = (s_LJ/d)**6.0_xp
!                 de = de + (q**2.0_xp - q)
!             end if
!         end do
!     end do

!     ! interactions between i and j images
!     !$omp parallel do private(j,d,q) reduction(+:de)
!     do j=1,npart*26
!         d = norm2(xi_new - x_PBC(j,:))
!         if (d .lt. rcut) then
!             q = (s_LJ/d)**6.0_xp
!             de = de + (q**2.0_xp - q)
!         end if
!     end do

!     ! finish
!     de = 4.0_dp*e0_LJ*de
!     end function erg_diff
