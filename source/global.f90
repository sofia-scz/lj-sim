module global
implicit none
integer, parameter :: dp = selected_real_kind(15, 307), &
                      xp = selected_real_kind(18, 4931), &
                      qp = selected_real_kind(33, 4931), &
                      int64 = selected_int_kind(18)
real(dp), parameter :: zero = 0.0d0, &
                       kB = 8.6173332620d-2, & 
                       amuangs_to_kglt = 1.660539066d0, &
                       mevangs_to_bar = 1.6021766208d3
integer :: npart
real(dp) :: mass, e0, sigma, rcut, xvar, vvar

! notes:
! kB in meV K is 8.6173303e-2 (NIST)
!
! LJ params
! neon
! e0 = 3.63242 meV 
! s = 2.76125 angs
! mass = 20.180
! expected a ~ 4.46 angs
!
! conversion 1 amu = 1.660539066e-24 g (NIST)
! c in atomic units = 2.9979e3 angs / fs  (NIST)
! Nav = 6.02214076e23
! pressure    1 meV/angs^3 = 1602.1766208 bar (berkeley)
end module global