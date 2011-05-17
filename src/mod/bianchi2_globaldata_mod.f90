!------------------------------------------------------------------------------
! bianchi2_globaldata_mod -- BIANCHI2 library global data class
!
!! Stores global data required in bianchi2 data routines.  Ideally one would 
!! avoid global data passed in this way, however it is required to get data
!! into the routines taken by the numerical recipes functions used in
!! computing bianchi2 templates.
!
!! @author A. N. Lasenby & J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 March 2006
!
! Revisions:
!   March 2006 - Written by Jason McEwen 
!------------------------------------------------------------------------------

module bianchi2_globaldata_mod

  use s2_types_mod

  integer, parameter :: b2gd_ntmax = 101
  integer, parameter :: b2gd_nvars = 4

  integer :: b2gd_it
  integer :: b2gd_nt

  real :: b2gd_treal(b2gd_ntmax)
  real :: b2gd_xreal(b2gd_ntmax)
!   real, allocatable :: b2gd_treal(:)
!   real, allocatable :: b2gd_xreal(:)

  real(s2_dp) :: b2gd_Omega_matter
  real(s2_dp) :: b2gd_Omega_Lambda
  real(s2_dp) :: b2gd_alpha, b2gd_Bianchi_h
  real(s2_dp) :: b2gd_RH_start
  real(s2_dp) :: b2gd_ze
  real(s2_dp) :: b2gd_theta_0, b2gd_phi_0     ! ANL's angle convention.

  real(s2_dp) :: b2gd_xarr(b2gd_nvars, b2gd_ntmax)
  real(s2_dp) :: b2gd_tarr(b2gd_ntmax)
!   real(s2_dp), allocatable :: b2gd_tarr(:)
!   real(s2_dp), allocatable :: b2gd_xarr(:, :)

  real(s2_dp) :: b2gd_tstop
  real(s2_dp) :: b2gd_deltat

end module bianchi2_globaldata_mod
