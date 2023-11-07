module global
implicit none
integer, parameter :: dp = selected_real_kind(15, 307), &
                      xp = selected_real_kind(18, 4931), &
                      qp = selected_real_kind(33, 4931), &
                      int64 = selected_int_kind(18)
real(dp), parameter :: zero = 0.0_dp, &
                       kB = 8.6173332620e-2, & 
                       amuang_to_kglt = 1.660539066, &
                       evang_to_bar = 1.6021766208e3
integer :: nPart
real(dp) :: mass, e0, sigma, frcut
! notes:
! kB in eV K is 8.6173303e-5
!
! LJ params
! neon
! e0 = 3.63242 meV 
! s = 2.76125 angs
! mass = 20.180
! expected a ~ 4.46 angs
!
! conversion 1 amu = 1.073544e-6 meV = 1.660539066e-24 g
! c in atomic units = 2.9979e3 angs / fs
! Nav = 6.02214076e23
! pressure    1 meV/angs^3 = 1602.1766208 bar
end module global