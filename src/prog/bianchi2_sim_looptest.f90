program bianchi2_sim_looptest

  use s2_types_mod
  use bianchi2_sky_mod

  implicit none

  type(bianchi2_sky) :: b
  integer :: i

  real(s2_dp) :: omega_matter, omega_lambda, h, zE, wH
  integer :: nside
  logical :: rhand

  ! Set parameters.
  omega_matter = 0.3d0
  omega_lambda = 0.6d0
  h = 1d-2
  zE = 1d3
  wH = 1d-10
  nside = 128
  rhand = .true.

  do i = 1,10000

    ! Initialise bianchi2 object.
    b = bianchi2_sky_init(omega_matter, omega_lambda, h, zE, wH, rhand, &
      nside)
  
    ! Free bianchi2 object.
    call bianchi2_sky_free(b)

  end do

end program bianchi2_sim_looptest