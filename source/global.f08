module global
use iso_fortran_env, only:  dp => real64, &
                            int64
implicit none
real(dp), parameter ::      pi = acos(-1.0d0), &
                            kB = 8.6173332620d-5, & 
                            amuangs_to_kglt = 1.660539066d0, &
                            evangs_to_bar = 1.6021766208d6, &
                            amu_speed = 1.036426965271228d2
integer :: npart
real(dp) :: mass, e0, sigma, rcut, xvar, vvar

! notes:
! 
! kB in eV K is 8.6173303e-5 (NIST)
! conversion 1 amu = 1.660539066e-24 g (NIST)
! c in atomic units = 2.9979e3 angs / fs  (NIST)
! Nav = 6.02214076e23
! pressure    1 meV/angs^3 = 1602.1766208 bar (berkeley)
end module global