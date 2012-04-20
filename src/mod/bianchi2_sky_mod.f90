!------------------------------------------------------------------------------
! bianchi2_sky_mod -- BIANCHI2 library sky class
!
!! Provides functionality to simulate a Bianchi2 VII_h model of the CMB, 
!! incorporating a cosmological constant.
!! Uses the s2_sky module to create healpix sky maps.
!
!! @author A. N. Lasenby & J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 October 2005
!
! Revisions:
!   October 2005 - Written by Anthony Lasenby & Jason McEwen 
!                  ANL has written all numerical code to compute Bianchi
!                  induced temperature fluctuations.  JDM has merely compiled
!                  this into a fortran 90 package that interfaces with 
!                  Healpix using the s2 library.
!------------------------------------------------------------------------------

module bianchi2_sky_mod

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod
  use bianchi2_error_mod

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    bianchi2_sky_init, &
    bianchi2_sky_free, &
    bianchi2_sky_write, &
    bianchi2_sky_rotate, &
    bianchi2_sky_param_write, &
    bianchi2_sky_compute_map, &
    bianchi2_sky_compute_alm, &
    bianchi2_sky_get_sky, &
    bianchi2_sky_get_alm, &
    bianchi2_sky_apply_beam


  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  ! None.


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  ! CMB T to convert Delta_T / T map to just Delta_T map.
#ifdef MILLIK
  ! Produce maps in units of mK.
  real(s2_dp), parameter :: BIANCHI2_CMB_T = 2.725d3
#else
  ! Produce Delta_T / T maps.
  real(s2_dp), parameter :: BIANCHI2_CMB_T = 1d0 
#endif


  !---------------------------------------
  ! Data types
  !---------------------------------------

  !! - init: Initialisation status.
  !! - omega_matter: Matter density (input parameter).
  !! - omega_lambda: Lambda density (input parameter).
  !! - h: Bianchi h parameter that controls `spiralness' (related to 
  !!   characteristic wavelength over which basis vectors change 
  !!   orientation) (input parameter).
  !! - zE: Red shift (input parameter).
  !! - wH: Normalised vorticity (normalised to Hubble constant) (input 
  !!   parameter).
  !! - rhand: Logical specifying handedness of map (true=right).
  !! - sky: Simulated map.

  type, public :: bianchi2_sky
     private
     logical :: init = .false.
     real(s2_dp) :: omega_matter = 0.0d0
     real(s2_dp) :: omega_lambda = 0.0d0
     real(s2_dp) :: h = 0.0d0
     real(s2_dp) :: zE = 0.0d0
     real(s2_dp) :: wH = 0.0d0
     logical :: rhand = .true.
     type(s2_sky) :: sky
  end type bianchi2_sky


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! bianchi2_sky_init
    !
    !! Initialise bianchi2 object by performing a bianchi2 simulation that
    !! incorporates a cosmological constant.
    !!
    !! Variables:
    !!  - omega_matter_in: Input omega_matter parameter (see bianchi2 data
    !!    type for explanation).
    !!  - omega_lambda_in: Input omega_lambda parameter (see bianchi2 data
    !!    type for explanation).
    !!  - h: Input h parameter (see bianchi2 data type for explanation).
    !!  - zE_in: Input zE parameter (see bianchi2 data type for explanation).
    !!  - wH: Input wH parameter (see bianchi2 data type for explanation).
    !!  - rhand: Logical to specify handedness of map.
    !!  - nside: Nside of Healpix map to generate.
    !!  - b: Initialised bianchi2 object with simulated map calculated.
    !
    !! @author A.N. Lasenby and J. D. McEwen
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Bianci computation part written by ANL.  Compiled into
    !                  init function and interfaced with s2 library by JDM.
    !--------------------------------------------------------------------------

    function bianchi2_sky_init(omega_matter_in, omega_lambda_in, h, zE_in, wH, rhand, &
         nside) result(b)

      use bianchi2_globaldata_mod
      use pix_tools, only: pix2ang_ring, nside2npix, in_ring 

      real(s2_dp), intent(in) :: omega_matter_in, omega_lambda_in, h, zE_in, wH
      logical, intent(in) :: rhand
      integer, intent(in) :: nside
      type(bianchi2_sky) :: b

      real(s2_dp) :: thetaOB, phiOB     ! Observing angles.
      real(s2_dp) :: theta0, phi0       ! Barrow's angle convention. 
  
      integer :: nt_use
      real(s2_dp) :: tstop_use
      real(s2_dp) :: tau_needed
      real(s2_dp) :: fact, U10_req, Phi1
      real(s2_dp) :: R_final, RH_final, cos_bit_final, sin_bit_final, other_bit
      real(s2_dp) :: Lambda, final_dens, t0

      integer :: npix, ipix, fail
      real(s2_dp), allocatable :: map(:)
      real(s2_sp) :: handedness_sign = +1d0

      real(s2_dp), parameter :: PHI_CENTRE = 0.0d0  ! Extract extire rings.
      real(s2_dp), parameter :: PHI_DIFF = PI
      integer :: iring, nring, max_pix_per_ring, npixring
      integer, allocatable :: ipixring(:)

      ! Check object not already initialised.
      if(b%init) then
        call bianchi2_error(BIANCHI2_ERROR_INIT, 'bianchi2_sky_init')
        return
      end if

      ! Initialise parameters passed as arguments.
      b%omega_matter = omega_matter_in
      b%omega_lambda = omega_lambda_in
      b%h = h
      b%zE = zE_in
      b%wH = wH
      b%rhand = rhand

      ! Initialise healpix map settings.
      npix = nside2npix(nside)
      allocate(map(0:npix-1), stat=fail)
      map = 0d0
      if(fail /= 0) then
         call bianchi2_error(BIANCHI2_ERROR_MEM_ALLOC_FAIL, &
           'bianchi2_sky_init')
      end if

      ! Set handedness sign.
      if(b%rhand) then
         handedness_sign = 1d0
      else
         handedness_sign = -1d0
      end if

!       ! Allocate global data.
!       allocate(b2gd_treal(1:b2gd_ntmax), stat=fail)
!       allocate(b2gd_xreal(1:b2gd_ntmax), stat=fail)
!       allocate(b2gd_tarr(1:b2gd_ntmax), stat=fail)
!       allocate(b2gd_xarr(1:b2gd_nvars, 1:b2gd_ntmax), stat=fail)
!       if(fail /= 0) then
!          call bianchi2_error(BIANCHI2_ERROR_MEM_ALLOC_FAIL, &
!            'bianchi2_sky_init')
!       end if

      ! Set global data.
      b2gd_Bianchi_h = b%h
      b2gd_Omega_matter = b%omega_matter
      b2gd_Omega_Lambda = b%omega_lambda
      b2gd_ze = b%zE
      b2gd_alpha = sqrt(b%h)
      b2gd_RH_start = sqrt(-b2gd_alpha**2/(b%omega_matter+b%omega_lambda-1d0))

      ! Initialise other variables.
      nt_use = 3
      b2gd_alpha = sqrt(b%h)
      b2gd_RH_start = sqrt(-b2gd_alpha**2/(b%omega_matter+b%omega_lambda-1d0))
      call get_tau(tau_needed)
      tstop_use=-tau_needed
      Lambda = 3*b2gd_RH_start**2*b%omega_lambda

      ! Set ring parameter values.
      nring = 4*nside - 1
      max_pix_per_ring = 4*nside

      ! Allocate space for ring pixel indices.
      allocate(ipixring(0:max_pix_per_ring-1), stat=fail)
      if(fail /= 0) then
         call bianchi2_error(BIANCHI2_ERROR_MEM_ALLOC_FAIL, &
           'bianchi2_sky_init')
      end if
      ipixring = 0

      ! Compute map value for each pixel in each ring.
      ! Do ring at a time since A and B functions of theta so once need
      ! to compute once for each ring.
      do iring = 1, nring    ! Checked this and ring is indeed indexed from 1.

!          if(mod(iring,nring/4) == 0) then
!             write(*,'(a,f5.1,a)') ' Percent complete: ', &
!                  iring/real(nring,s2_sp)*100e0, '%'
!          end if

         ! Get ring_pix in order to determine theta of ring.
         call in_ring(nside, iring, PHI_CENTRE, PHI_DIFF, ipixring, &
              & npixring, nest=0) ! Ring scheme!!!!!!!!!

         ! Calculate thetaOB for current ring.
         call pix2ang_ring(nside, ipixring(0), thetaOB, phiOB)

         ! Set photon theta angle from observation angle.
         theta0 = PI - thetaOB

         ! Set global data.
         ! (Compute theta_0 for ANL's convention.)
         b2gd_theta_0 = theta0 - pi/2d0

         ! Compute terms for constant theta.
         call get_results(tstop_use,R_final,RH_final,cos_bit_final,sin_bit_final)
         fact = sqrt(9.D0*b2gd_alpha**2+1.D0)*sqrt(1.D0+b2gd_alpha**2)* &
                sqrt(1.D0-b%omega_matter-b%omega_lambda)**3.D0/b2gd_alpha**3 &
                /b%omega_matter/6.D0
         U10_req = b%wH/fact

         ! Compute map valus for all pixels in current ring.
         do ipix = 0,npixring-1

            ! Calculate thetaOB for current ring.
            call pix2ang_ring(nside, ipixring(ipix), thetaOB, phiOB)

            ! Set photon phi angle from observation angle.
            phi0 = phiOB - PI

            ! Account for handedness 
            ! (also component outside loop to account for handedness).
            phi0 = phi0 * handedness_sign

            !  Compute phi_0 for ANL's convention.
            b2gd_phi_0 = phi0 + 3d0/4d0*pi

            ! Get first contribution.
            Phi1 = U10_req/R_final*(sin_bit_final*sin(b2gd_phi_0)+cos_bit_final*cos(b2gd_phi_0))
            map(ipixring(ipix)) = -Phi1/(1+b2gd_ze)*2d0   ! JDM: Need factor 2 to make agree

            ! Now get end contribution.
            final_dens = (3d0*RH_final**2-Lambda*R_final**2-3d0*b2gd_alpha**2)/(R_final**2)
            t0 = 3.D0/final_dens/R_final**5*cos((b2gd_phi_0*b2gd_alpha-b2gd_alpha*tstop_use- &
                dlog(cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0+exp(-2.D0*b2gd_alpha*tstop_use)- &
                exp(-2.D0*b2gd_alpha*tstop_use)*cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0))/b2gd_alpha)* &
                exp(b2gd_alpha*tstop_use)*cos(b2gd_theta_0)/(-1.D0-sin(b2gd_theta_0)+ &
                exp(2.D0*b2gd_alpha*tstop_use)*sin(b2gd_theta_0)-exp(2.D0*b2gd_alpha*tstop_use))*b2gd_alpha+ &
                1.D0/final_dens/R_final**5*sin((b2gd_phi_0*b2gd_alpha-b2gd_alpha*tstop_use- &
                dlog(cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0+exp(-2.D0*b2gd_alpha*tstop_use)- & 
                exp(-2.D0*b2gd_alpha*tstop_use)*cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0))/b2gd_alpha)* &
                exp(b2gd_alpha*tstop_use)*cos(b2gd_theta_0)/(-1.D0-sin(b2gd_theta_0)+ &
                exp(2.D0*b2gd_alpha*tstop_use)*sin(b2gd_theta_0)-exp(2.D0*b2gd_alpha*tstop_use))
            other_bit = U10_req*t0

            map(ipixring(ipix)) = map(ipixring(ipix)) - other_bit/(1+b2gd_ze)*2d0 

         end do

      end do

!       write(*,'(a)') ' Percent complete: 100.0%'

      ! Account for handedness.
      map = handedness_sign * map

      ! Convert Delta_T/T map computed to Delta_T map.
      map = BIANCHI2_CMB_T * map

      ! Initialise sky object with map.
      b%sky = s2_sky_init(real(map,s2_sp), nside, S2_SKY_RING)

      ! Set initialised status.
      b%init = .true.

!       ! Free global data.
!       deallocate(b2gd_treal)
!       deallocate(b2gd_xreal)
!       deallocate(b2gd_tarr)
!       deallocate(b2gd_xarr)

      ! Free memory.
      deallocate(map)
      deallocate(ipixring)

    end function bianchi2_sky_init



    !--------------------------------------------------------------------------
    ! bianchi2_sky_free
    !
    !! Free all memory associated with a bianchi2 object.
    !!
    !! Variables:
    !!  - b: Bianchi2 object to free.
    !
    !! @author J. D. McEwen
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi2_sky_free(b)

      type(bianchi2_sky), intent(inout) :: b

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi2_error(BIANCHI2_ERROR_NOT_INIT, 'bianchi2_sky_free')
      end if 

      ! Free sky.
      call s2_sky_free(b%sky)

      ! Reset attributes.
      b%omega_matter = 0.0d0
      b%omega_lambda = 0.0d0
      b%h = 0.0d0
      b%zE = 0.0d0
      b%wH = 0.0d0
      b%rhand = .true.
      b%init = .false.

    end subroutine bianchi2_sky_free


    !--------------------------------------------------------------------------
    ! bianchi2_sky_param_write
    !
    !! Write parameters of the Bianchi2 simulation to standard output.
    !!
    !! Variables:
    !!  - b: Bianchi2 object containing parameter attributes to write.
    !
    !! @author J. D. McEwen
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi2_sky_param_write(b)

      type(bianchi2_sky), intent(in) :: b

      write(*,'(a,e12.5)') ' b%omega_matter: ', b%omega_matter
      write(*,'(a,e12.5)') ' b%omega_lmabda: ', b%omega_lambda
      write(*,'(a,e23.5)') ' b%h: ', b%h
      write(*,'(a,e22.5)') ' b%zE: ', b%zE
      write(*,'(a,e22.5)') ' b%wH: ', b%wH
      write(*,'(a,l19)') ' b%rhand: ', b%rhand

    end subroutine bianchi2_sky_param_write


    !--------------------------------------------------------------------------
    ! bianchi2_sky_compute_map
    !
    !! Compute the map of the bianchi2 sky, assuming the alms are already
    !! defined.
    !!
    !! Variables:
    !!  - b: Bianchi2 object containing alms to compute map of.
    !!  - nside: Healpix nside to compute map at.
    !
    !! @author J. D. McEwen
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi2_sky_compute_map(b, nside)

      type(bianchi2_sky), intent(inout) :: b
      integer, intent(in) :: nside

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi2_error(BIANCHI2_ERROR_NOT_INIT, &
          'bianchi2_sky_compute_map')
      end if

      ! Compute bianchi2 sky alms.
      call s2_sky_compute_map(b%sky, nside)

    end subroutine bianchi2_sky_compute_map


    !--------------------------------------------------------------------------
    ! bianchi2_sky_compute_alm
    !
    !! Compute the alm of the bianchi2 sky, assuming the map is already
    !! defined.
    !!
    !! Variables:
    !!  - b: Bianchi2 object containing sky to compute alms of.
    !!  - lmax: Maximum harmonic l to consider when computing alms.
    !!  - mmax: Maximum harmonic m to consider when computing alms.
    !
    !! @author J. D. McEwen
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi2_sky_compute_alm(b, lmax, mmax)

      type(bianchi2_sky), intent(inout) :: b
      integer, intent(in) :: lmax, mmax

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi2_error(BIANCHI2_ERROR_NOT_INIT, &
          'bianchi2_sky_compute_alm')
      end if

      ! Compute bianchi2 sky alms.
      call s2_sky_compute_alm(b%sky, lmax, mmax)

    end subroutine bianchi2_sky_compute_alm


    !--------------------------------------------------------------------------
    ! bianchi2_sky_get_sky
    !
    !! Get sky variable from the passed bianchi2 object.
    !!
    !! Notes:
    !!   - Initialises a new sky as a copy of the bianchi2 sky.
    !!   - The returned sky is subsequently independed of the sky stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - b: Bianchi2 object containing sky to get.
    !!   - sky: Object sky variable returned.
    !
    !! @author J. D. McEwen 
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function bianchi2_sky_get_sky(b) result(sky)

      type(bianchi2_sky), intent(in) :: b
      type(s2_sky) :: sky

      ! Check object initialised.
      if(.not. b%init) then
        call s2_error(BIANCHI2_ERROR_NOT_INIT, 'bianchi2_sky_get_sky')
      end if 

      ! Make a copy for the returned sky.
      sky = s2_sky_init(b%sky)

    end function bianchi2_sky_get_sky

    
    !--------------------------------------------------------------------------
    ! bianchi2_sky_get_alm
    !
    !! Get alms contained in a bianchi2 sky object.  Alms must already be
    !! computed.
    !!
    !! Notes:
    !!   - Error occurs if alms are not already computed.
    !!
    !! Variables:
    !!   - b: Bianchi2 object containing sky that in turn contained the alms
    !!     to get.
    !!   - alm: Alms of simulated bianchi map extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi2_sky_get_alm(b, alm)

      type(bianchi2_sky), intent(in) :: b
      complex(s2_spc), intent(out) :: alm(:,:)

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi2_error(BIANCHI2_ERROR_NOT_INIT, 'bianchi2_sky_get_alm')
      end if 

      ! Get alms from sky.
      call s2_sky_get_alm(b%sky, alm)

    end subroutine bianchi2_sky_get_alm


    !--------------------------------------------------------------------------
    ! bianchi2_sky_apply_beam
    !
    !! Apply Gaussian beam with specified FWHM.  (FWHM must be passed in
    !! arcmin.)
    !!
    !! Notes:
    !!   - Error occurs if alms are not already computed.
    !!
    !! Variables:
    !!   - b: Bianchi2 object containing sky that is to be comvolved with
    !!     the beam.
    !!   - fwhm: Gaussian beam FWHM to use (specified in arcmin).
    !!   - lmax: Maximum harmonic l to consider.
    !
    !! @author J. D. McEwen
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi2_sky_apply_beam(b, fwhm, lmax)

      use s2_pl_mod
      use s2_vect_mod, only: s2_vect_arcmin_to_rad

      type(bianchi2_sky), intent(inout) :: b
      real(s2_sp), intent(in) :: fwhm
      integer, intent(in) :: lmax

      real(s2_sp) :: fwhm_rad
      type(s2_pl) :: beam

      ! Convert beam_fwhm to radians.
      fwhm_rad = s2_vect_arcmin_to_rad(fwhm)

      ! Create beam.
      beam = s2_pl_init_guassian(fwhm_rad, lmax)

      ! Apply beam.
      call s2_sky_conv(b%sky, beam)

      ! Free beam.
      call s2_pl_free(beam)

    end subroutine bianchi2_sky_apply_beam


    !--------------------------------------------------------------------------
    ! bianchi2_sky_write
    !
    !! Write the bianchi2 simulated sky to a file. 
    !!
    !! Variables:
    !!  - b: Bianchi2 object to save sky of.
    !!  - filename: Name of output file.
    !!  - file_type: Type of output file, either fits map or sky file 
    !!    (see s2_sky_mod for more details).  Integer flag that may be 
    !!    either S2_SKY_FILE_TYPE_MAP or S2_SKY_FILE_TYPE_SKY.
    !!  - [comment]: Optional comment to append to file header.
    !
    !! @author J. D. McEwen
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine bianchi2_sky_write(b, filename, file_type, comment)

      type(bianchi2_sky), intent(inout) :: b
      character(len=*), intent(in) :: filename
      integer, intent(in) :: file_type
      character(len=*), intent(in), optional :: comment

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi2_error(BIANCHI2_ERROR_NOT_INIT, 'bianchi2_sky_write')
      end if 

      select case(file_type)

         case(S2_SKY_FILE_TYPE_MAP)
            call s2_sky_write_map_file(b%sky, filename, comment)

         case(S2_SKY_FILE_TYPE_SKY)
            call s2_sky_io_fits_write(filename, b%sky, comment)

         case default
            call s2_error(S2_ERROR_SKY_FILE_INVALID, 'bianchi2_sky_write', &
              comment_add='Invalid file type specifier')

      end select

    end subroutine bianchi2_sky_write


    !--------------------------------------------------------------------------
    ! bianchi2_sky_rotate
    !
    !! Rotate the simulated sky of the bianchi2 object.  
    !!
    !! Variables:
    !!  - b: Bianchi2 object containing sky to be rotated.
    !!  - alpha: Alpha Euler angle of the rotation.
    !!  - beta: Beta Euler angle of the rotation.
    !!  - gamma: Gamma Euler angle of the rotation.
    !!  - lmax: Maximum harmonic l to consider.
    !!  - nside: Healpix nside to compute map at.
    !!  - rotation_alm : put .false. in bianchi2_sim.f90 if
    !     one wants to perform the rotation pixel by pixel
    !
    !! @author J. D. McEwen
    !! @version 0.1 October 2005
    !
    ! Revisions:
    !   October 2005 - Written by Jason McEwen
    !   April 2012 - Modifications by Thibaut Josset
    !                Added option to perform rotation in harmonic space.
    !--------------------------------------------------------------------------

    subroutine bianchi2_sky_rotate(b, alpha, beta, gamma, lmax, nside,&
 rotation_alm)

      use s2_dl_mod, only: s2_dl_beta_operator

      type(bianchi2_sky), intent(inout) :: b
      real(s2_sp), intent(in), optional :: alpha, beta, gamma

      logical, intent(in) :: rotation_alm
      integer, intent(in) :: lmax, nside
      integer :: l, m
      integer :: fail=10
      complex(s2_spc), allocatable :: alm(:,:)
      complex(s2_spc), allocatable :: alm_rotate(:,:)
      real(s2_dp), pointer :: dl(:,:) => null()
      complex(s2_dpc) :: Dm_p1, Dm_m1
      real(s2_sp), parameter :: ZERO_TOL = 1d-4
      complex(s2_dpc) :: icmpx
      !set icmpx=sqrt(-1)
      icmpx = cmplx(0d0,1d0)

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi2_error(BIANCHI2_ERROR_NOT_INIT, 'bianchi2_sky_rotate')
      end if
      

      ! Rotate if Euler angles present and at least one angle non-zero
      if(present(alpha) .and. present(beta) .and. present(gamma) &
           .and. (abs(alpha)+abs(beta)+abs(gamma) > ZERO_TOL) ) then
         
         !Choose between rotation pixel by pixel or rotation of the alm
         if (rotation_alm == .false.) then

            write(*,*) 'Rotation pixel by pixel with s2_sky_rotate'
            call s2_sky_rotate(b%sky, alpha, beta, gamma)

         else

            write(*,*) 'Rotation of the alm'
            allocate(dl(-lmax:lmax,-lmax:lmax),stat=fail)
            allocate(alm(0:lmax,0:lmax), stat=fail)
            allocate(alm_rotate(0:lmax,0:lmax), stat=fail)
            if(fail/=0) then
               call bianchi2_error(BIANCHI2_ERROR_MEM_ALLOC_FAIL, &
                    'bianchi_sky_rotate')
            endif
            alm_rotate = 0e0

            !get alm
            call bianchi2_sky_compute_alm(b,lmax,lmax)
            call bianchi2_sky_get_alm(b, alm)

            ! perform rotation in harmonic space,&
            ! noting Bianchi alms only non zero for m=+-1)
            do l = 1,lmax

               call s2_dl_beta_operator(dl, real(beta,s2_dp), l)
            
               do m = 0,l


                  !Calculation of  Dl,m,+-1
                  Dm_p1 = exp(-icmpx*m*alpha) * dl(m, 1) * exp(-icmpx*gamma)
                  Dm_m1 = exp(-icmpx*m*alpha) * dl(m,-1) * exp( icmpx*gamma)

                  !Rotation of the alm
                  alm_rotate(l,m) = - Dm_m1*conjg(alm(l,1)) + Dm_p1*alm(l,1)
             
               end do

           end do
           
           !Make free b%sky
           call s2_sky_free(b%sky)

           ! Compute sky object with rotated alms
           b%sky = s2_sky_init(alm_rotate, lmax, lmax, nside)
           
           call bianchi2_sky_compute_map(b,nside)

           ! Dellocate memory used for rotating alms
           deallocate(alm_rotate)
           deallocate(alm)
           deallocate(dl)

         end if

       endif

    end subroutine bianchi2_sky_rotate


    !--------------------------------------------------------------------------
    ! Routines to compute Bianchi induced fluctuations
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! get_results
    !
    !! This returns the values found at the end of the integral
    !
    !! @author A. N. Lasenby
    !! @version 08/10/05
    !
    ! Revisions:
    !   October 2005 - Incorporated into library by JDM
    !--------------------------------------------------------------------------

    subroutine get_results(tstop_use, R_final, RH_final, &
         cos_bit_final, sin_bit_final)

      use s2_types_mod
      use bianchi2_globaldata_mod

      implicit none

      real(s2_dp), intent(in) :: tstop_use
      real(s2_dp), intent(out) :: R_final,RH_final,cos_bit_final,sin_bit_final

      integer, parameter :: neq = 4

      real(s2_dp) :: vals(b2gd_nvars)
      real(s2_dp) :: warr(21*neq+28)
      real(s2_dp) :: tstart,tol
      real(s2_dp) :: tmin,tmax,t,tend    
      real(s2_dp) :: R,RH
#ifdef NAGI8
      integer(kind=8) :: neqn,irelab,ifail
#else
      integer :: neqn,irelab,ifail
#endif

      character(len=1) :: relabs
      external d02cjw

      b2gd_tstop=tstop_use
      tstart=0
	
      b2gd_nt=3

      b2gd_deltat=(b2gd_tstop-tstart)/(b2gd_nt-1)

      tmin=tstart
      tmax=b2gd_tstop

      tol=1d-12

      t=tstart

      ! Next bit sets up initial value for R

      vals(1)=1

      !	Next bit sets up initial value for RH

      vals(2)=b2gd_RH_start

      vals(3)=0d0
      vals(4)=0d0

      t=tstart
      tend=b2gd_tstop
      neqn=neq
      irelab=0
      ifail=1
      b2gd_it=b2gd_nt-2
      relabs='D'

      call d02cjf(t,tend,neqn,vals,fcn,tol,relabs,out,d02cjw,warr,ifail)

      ! Check success.
      if((ifail/=0).and.(ifail/=1)) then
        call bianchi2_error(BIANCHI2_ERROR_SKY_NUM_FAIL, 'get_results', &
          comment_add='NAG routine d02cjf failed')
      end if
      
      do b2gd_it=1,b2gd_nt
         b2gd_treal(b2gd_it)=b2gd_tarr(b2gd_it)
         R=b2gd_xarr(1,b2gd_it)
         RH=b2gd_xarr(2,b2gd_it)
         t=b2gd_tarr(b2gd_it)
         b2gd_xreal(b2gd_it)=R
      end do

      R_final=b2gd_xarr(1,b2gd_nt)
      RH_final=b2gd_xarr(2,b2gd_nt)	
      cos_bit_final=b2gd_xarr(3,b2gd_nt)	
      sin_bit_final=b2gd_xarr(4,b2gd_nt)

    end subroutine get_results


    !--------------------------------------------------------------------------
    ! fcn
    !
    !! @author A. N. Lasenby
    !! @version 08/10/05
    !
    ! Revisions:
    !   October 2005 - Incorporated into library by JDM
    !--------------------------------------------------------------------------

    subroutine fcn(t,vars,ff)

      use s2_types_mod
      use bianchi2_globaldata_mod

      implicit none
      integer, parameter :: neq = 4
      integer, parameter :: nvars = 4

      real(s2_dp), intent(in) :: t
      real(s2_dp), intent(in) :: vars(nvars)
      real(s2_dp), intent(out) :: ff(nvars)

      real(s2_dp) ::tau
      real(s2_dp) :: R,RH,Lambda
      real(s2_dp) :: arg,fact

      R=vars(1)
      RH=vars(2)
      tau=t

      Lambda=3*b2gd_RH_start**2*b2gd_Omega_Lambda

      ff(1)=R*RH

      ff(2)=-RH**2.D0/2.D0+Lambda*R**2.D0/2.D0+b2gd_alpha**2/2.D0

      arg=(b2gd_alpha*tau+dlog(cos(0.3141592653589793D1/4.D0+ &
           b2gd_theta_0/2.D0)**2.D0+exp(-2.D0*b2gd_alpha*tau)-&
           exp(-2.D0*b2gd_alpha*tau)*cos(0.3141592653589793D1/4.D0+ &
           b2gd_theta_0/2.D0)**2.D0))/b2gd_alpha

      fact=-2.D0/R**2*exp(b2gd_alpha*tau)*cos(b2gd_theta_0)/(-1.D0-sin(b2gd_theta_0)+ &
           exp(2.D0*b2gd_alpha*tau)*sin(b2gd_theta_0)-exp(2.D0*b2gd_alpha*tau))**2.D0* &
           (1.D0+sin(b2gd_theta_0)+exp(2.D0*b2gd_alpha*tau)*sin(b2gd_theta_0)- &
           exp(2.D0*b2gd_alpha*tau))

      ff(3)=fact*cos(arg)
      
      ff(4)=fact*sin(arg)

    end subroutine fcn


    !--------------------------------------------------------------------------
    ! out
    !
    !! @author A. N. Lasenby
    !! @version 08/10/05
    !
    ! Revisions:
    !   October 2005 - Incorporated into library by JDM
    !--------------------------------------------------------------------------

    subroutine out(t,vals)

      use s2_types_mod
      use bianchi2_globaldata_mod

      implicit none

      real(s2_dp), intent(inout) :: t
      real(s2_dp), intent(in) :: vals(b2gd_nvars)

      real(s2_dp) :: deriv(b2gd_nvars)
      integer :: iel, ival

      iel=b2gd_nt-(b2gd_it+1)
      b2gd_tarr(iel)=t
      do ival=1,b2gd_nvars
         b2gd_xarr(ival,iel)=vals(ival)
      end do
      
      call fcn(t,vals,deriv)

      t=b2gd_tstop-b2gd_it*b2gd_deltat
      b2gd_it=b2gd_it-1

    end subroutine out


    !--------------------------------------------------------------------------
    ! get_tau
    !
    !! @author A. N. Lasenby
    !! @version 08/10/05
    !
    ! Revisions:
    !   October 2005 - Incorporated into library by JDM
    !--------------------------------------------------------------------------

    subroutine get_tau(tau_needed)

      use s2_types_mod
      use bianchi2_globaldata_mod

      implicit none

      real(s2_dp), intent(out) :: tau_needed

      real(s2_dp) :: aa, ABSERR, bb, EPSABS, EPSREL, RESULT_VAL

#ifdef NAGI8
      integer(kind=8) :: IFAIL
      integer(kind=8), parameter :: LIWORK = 10000
      integer(kind=8), parameter :: LWORK = 10000
      integer(kind=8) :: IWORK(LIWORK)
#else
      integer :: IFAIL
      integer, parameter :: LIWORK = 10000
      integer, parameter :: LWORK = 10000
      integer :: IWORK(LIWORK)
#endif

      real(s2_dp) :: WORK(LWORK)

      ifail=0
      epsabs=1d-6
      epsrel=1d-8

      aa=1d0/sqrt(1+b2gd_ze)
      bb=1d0
      call d01ajf(F1r,aa,bb,EPSABS,EPSREL,RESULT_VAL,ABSERR,WORK,LWORK,&
        IWORK,LIWORK,IFAIL)
	
      tau_needed=result_val

    end subroutine get_tau


    !--------------------------------------------------------------------------
    ! F1r
    !
    !! @author A. N. Lasenby
    !! @version 08/10/05
    !
    ! Revisions:
    !   October 2005 - Incorporated into library by JDM
    !--------------------------------------------------------------------------

    function F1r(a)

      use s2_types_mod
      use bianchi2_globaldata_mod

      implicit none

      real(s2_dp) :: F1r
      real(s2_dp), intent(in) :: a

      real(s2_dp) :: t0

      t0 = -b2gd_alpha**2*(b2gd_Omega_Lambda*a**6-a**2*b2gd_Omega_matter-a**2* &
           b2gd_Omega_Lambda+a**2+b2gd_Omega_matter)/(b2gd_Omega_matter+b2gd_Omega_Lambda-1)

      F1r=2d0/sqrt(t0)

    end function F1r


end module bianchi2_sky_mod
