!==============================================================================
!> \mainpage Bianchi2
!!
!! The Bianchi2 package provides functionality to simulate Bianchi
!! Type VIIh induced temperature fluctuations in CMB maps of a
!! universe with shear and rotation. The implementation is based on
!! the solutions to the Bianchi models derived by Anthony Lasenby (not
!! yet published) that incorporate a cosmological
!! constant. Functionality is provided to compute the induced
!! fluctuations on the sphere directly in either real or harmonic
!! space.
!!
!! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!! \authors <a href="http://www.google.co.uk/search?q=Anthony%20Lasenby">
!!          Anthony Lasenby</a>
!! \authors Thibaut Josset
!==============================================================================


!------------------------------------------------------------------------------
! bianchi2_sky_mod -- BIANCHI2 library sky class
!
!> Provides functionality to simulate a Bianchi2 VII_h model of the CMB, 
!! incorporating a cosmological constant.
!! Uses the s2_sky module to create healpix sky maps.
!!
!! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!! \authors <a href="http://www.google.co.uk/search?q=Anthony%20Lasenby">
!!          Anthony Lasenby</a>
!------------------------------------------------------------------------------

module bianchi2_sky_mod

  use s2_types_mod
  use s2_sky_mod
  use s2_error_mod
  use bianchi2_error_mod
  use omp_lib

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    bianchi2_sky_init, &
    bianchi2_sky_init_alm, &
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
  ! Logical to disable the lookup table.
  logical, parameter :: BIANCHI2_DISABLE_PLM1TABLE = .false.

  !---------------------------------------
  ! Data types
  !---------------------------------------

  ! Bianchi 2 sky object that contains parameters and simulated sky.
  type, public :: bianchi2_sky
     private
     ! Initialisation status.
     logical :: init = .false.
     ! Matter density (input parameter).
     real(s2_dp) :: omega_matter = 0.0d0
     ! Lambda density (input parameter).
     real(s2_dp) :: omega_lambda = 0.0d0
     ! Bianchi h parameter that controls `spiralness' (related to 
     ! characteristic wavelength over which basis vectors change 
     ! orientation) (input parameter).
     real(s2_dp) :: h = 0.0d0
     ! Red shift (input parameter).
     real(s2_dp) :: zE = 0.0d0
     ! wH: Normalised vorticity (normalised to Hubble constant) (input 
     ! parameter).
     real(s2_dp) :: wH = 0.0d0
     ! Logical specifying handedness of map (true=right).
     logical :: rhand = .true.
     ! Simulated map.
     type(s2_sky) :: sky
  end type bianchi2_sky


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! bianchi2_sky_init
    !
    !> Initialise bianchi2 object by performing a bianchi2 simulation that
    !! incorporates a cosmological constant.
    !!
    !!  \param[in] omega_matter_in Input omega_matter parameter (see bianchi2 data
    !!    type for explanation).
    !!  \param[in] omega_lambda_in Input omega_lambda parameter (see bianchi2 data
    !!    type for explanation).
    !!  \param[in] h Input h parameter (see bianchi2 data type for explanation).
    !!  \param[in] zE_in Input zE parameter (see bianchi2 data type for explanation).
    !!  \param[in] wH Input wH parameter (see bianchi2 data type for explanation).
    !!  \param[in] rhand Logical to specify handedness of map.
    !!  \param[in] nside Nside of Healpix map to generate.
    !!  \retval b Initialised bianchi2 object with simulated map calculated.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
    !! \authors <a href="http://www.google.co.uk/search?q=Anthony%20Lasenby">
    !!          Anthony Lasenby</a>
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

      !> Initialise parameters passed as arguments.
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

      !> Set handedness sign.
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

      !> Set global data.
      b2gd_Bianchi_h = b%h
      b2gd_Omega_matter = b%omega_matter
      b2gd_Omega_Lambda = b%omega_lambda
      b2gd_ze = b%zE

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
      
      do iring = 1, nring    !> Checked this and ring is indeed indexed from 1.

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

            ! Compute phi_0 for ANL's convention.
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
    ! bianchi2_sky_init_alm
    !
    !> Initialise bianchi2 object by performing a bianchi2 simulation that
    !! incorporates a cosmological constant.
    !!
    !!  \param[in] omega_matter_in Input omega_matter parameter (see bianchi2 data
    !!    type for explanation).
    !!  \param[in] omega_lambda_in: Input omega_lambda parameter (see bianchi2 data
    !!    type for explanation).
    !!  \param[in] h: Input h parameter (see bianchi2 data type for explanation).
    !!  \param[in] zE_in: Input zE parameter (see bianchi2 data type for explanation).
    !!  \param[in] wH: Input wH parameter (see bianchi2 data type for explanation).
    !!  \param[in] rhand: Logical to specify handedness of map.
    !!  \param[in] nside: Nside of Healpix map to generate.
    !!  \param[in] lmax: Maximum harmonic l to consider.
    !!  \param[in] Nuse : Number of terms in use in the integrations IA and IB.
    !!  \param[in] alpha: Alpha Euler angle of the rotation.
    !!  \param[in] beta: Beta Euler angle of the rotation.
    !!  \param[in] gamma: Gamma Euler angle of the rotation.
    !!  \retval b Initialised bianchi2 object with simulated map calculated.
    !!    
    !!  \authors Thibaut Josset
    !--------------------------------------------------------------------------

    function bianchi2_sky_init_alm(omega_matter_in, omega_lambda_in, h, zE_in, wH, rhand, &
         nside, lmax, Nuse, alpha, beta, gamma) result(b)

      use bianchi2_globaldata_mod
      use s2_dl_mod, only : s2_dl_beta_operator

      real(s2_dp), intent(in) :: omega_matter_in, omega_lambda_in, h, zE_in, wH
      logical, intent(in) :: rhand
      integer, intent(in) :: nside, lmax, Nuse
      type(bianchi2_sky) :: b
  
      real(s2_dp) :: tstop_use
      real(s2_dp) :: tau_needed
      real(s2_dp) :: fact, U10_req
      real(s2_dp) :: R_final, RH_final, cos_bit_final, sin_bit_final
      real(s2_dp) :: Lambda, final_dens

      integer :: npix, ipix, fail

      ! Parameters for computing alm.
      complex(s2_spc), allocatable :: alm(:,:)
      real(s2_dp) :: handedness_sign = +1d0, lsign = 1d0
      integer :: l, itheta
      real(s2_dp), allocatable :: A_grid(:), B_grid(:)  
      real(s2_dp) :: theta_inc
      real(s2_dp) :: C_sin, C_cos
      real(s2_dp) :: IA, IB
      
      ! Parameters for the rotation.
      real(s2_sp), intent(in), optional :: alpha, beta, gamma
      real(s2_sp), parameter :: ZERO_TOL = 1d-4
      real(s2_dp), pointer :: dl(:,:) => null()
      integer :: m
      complex(s2_spc), allocatable :: alm_rotate(:,:)
      complex(s2_dpc) :: Dm_p1, Dm_m1
      complex(s2_dpc) :: icmpx

      ! Set icmpx=sqrt(-1).
      icmpx = cmplx(0d0,1d0)

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

      ! Initialise healpix alm settings.
      allocate(alm(0:lmax,0:1), stat=fail)
      alm = cmplx(0d0, 0d0)
      if(fail /= 0) then
        call bianchi2_error(BIANCHI2_ERROR_MEM_ALLOC_FAIL, &
          'bianchi2_sky_init_alm')
      end if

      ! Set handedness sign.
      if(b%rhand) then
         handedness_sign = 1d0
      else
         handedness_sign = -1d0
      end if

      ! Set global data.
      b2gd_Bianchi_h = b%h
      b2gd_Omega_matter = b%omega_matter
      b2gd_Omega_Lambda = b%omega_lambda
      b2gd_ze = b%zE

      ! Initialise other variables.
      b2gd_alpha = sqrt(b%h)
      b2gd_RH_start = sqrt(-b2gd_alpha**2/(b%omega_matter+b%omega_lambda-1d0))
      call get_tau(tau_needed)
      tstop_use=-tau_needed    

      ! Calculate A(theta) and B(theta) terms.
      allocate(A_grid(0:Nuse-1), stat=fail)
      allocate(B_grid(0:Nuse-1), stat=fail)

      if(fail /= 0) then
         call bianchi2_error(BIANCHI2_ERROR_MEM_ALLOC_FAIL, &
            'bianchi2_sky_init_alm')
      end if

      A_grid = 0d0
      B_grid = 0d0
      b2gd_theta_0 = -pi/2d0
      theta_inc = (pi - 0d0) / real(Nuse, s2_dp)

!This does not work :
!$OMP PARALLEL default(none), copyin(b2gd_it, b2gd_nt, b2_gd_treal,b2gd_xreal,b2gd_theta_0,b2gd_xarr,b2gd_tarrb2gd_RH_start,b2gd_deltat,b2gd_ze,b2gd_tstop,b2gd_alpha,b2gd_Omega_matter,b2gd_Omega_Lambda,b2gd_RH_start) &
!$OMP PRIVATE(itheta,R_final,RH_final,cos_bit_final,sin_bit_final,Lambda,final_dens,C_sin,C_cos,theta_inc,tstop_use,fact,U10_req,b) &
!$OMP SHARED(Nuse,A_grid,B_grid)


!$OMP DO
      do itheta = 0, Nuse-1

         b2gd_theta_0=-pi/2d0+itheta*theta_inc
        
         call get_results(tstop_use,R_final,RH_final,cos_bit_final,sin_bit_final)
         fact = sqrt(9.D0*b2gd_alpha**2+1.D0)*sqrt(1.D0+b2gd_alpha**2)* &
                sqrt(1.D0-b%omega_matter-b%omega_lambda)**3.D0/b2gd_alpha**3 &
                /b%omega_matter/6.D0
         U10_req = b%wH/fact
         Lambda = 3*b2gd_RH_start**2*b%omega_lambda
         final_dens = (3d0*RH_final**2-Lambda*R_final**2-3d0*b2gd_alpha**2)/(R_final**2)

        C_sin=-U10_req/(1+b2gd_ze)*2d0*( &
                1/R_final*(sin_bit_final*cos(3d0/4d0*pi)-cos_bit_final*sin(3d0/4d0*pi)) &
                - 3.D0/final_dens/R_final**5 &
                  * exp(b2gd_alpha*tstop_use)*cos(b2gd_theta_0) &
                  / (-1.D0-sin(b2gd_theta_0) &
                     + exp(2.D0*b2gd_alpha*tstop_use)*sin(b2gd_theta_0) &
                     - exp(2.D0*b2gd_alpha*tstop_use)) &
                  * b2gd_alpha &
                  * sin((b2gd_alpha*tstop_use- &
                          dlog(cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0 &
                              + exp(-2.D0*b2gd_alpha*tstop_use) &
                              - exp(-2.D0*b2gd_alpha*tstop_use)*cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0) &
                         )/b2gd_alpha+3d0/4d0*pi) &
               + 1.D0/final_dens/R_final**5 &
                 *  exp(b2gd_alpha*tstop_use)*cos(b2gd_theta_0) &
                 / (-1.D0-sin(b2gd_theta_0) &
                    + exp(2.D0*b2gd_alpha*tstop_use)*sin(b2gd_theta_0) &
                    -exp(2.D0*b2gd_alpha*tstop_use)) &
                 * cos((-b2gd_alpha*tstop_use- &
                        dlog(cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0 &
                            + exp(-2.D0*b2gd_alpha*tstop_use) &
                            - exp(-2.D0*b2gd_alpha*tstop_use)*cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0) &
                       )/b2gd_alpha+3d0/4d0*pi) &
               )

        C_cos=-U10_req/(1+b2gd_ze)*2d0*( &
                1/R_final*(sin_bit_final*sin(3d0/4d0*pi)+cos_bit_final*cos(3d0/4d0*pi)) &
                + 3.D0/final_dens/R_final**5 &
                  * exp(b2gd_alpha*tstop_use)*cos(b2gd_theta_0) &
                  / (-1.D0-sin(b2gd_theta_0) &
                     + exp(2.D0*b2gd_alpha*tstop_use)*sin(b2gd_theta_0) &
                     - exp(2.D0*b2gd_alpha*tstop_use)) &
                  * b2gd_alpha &
                  * cos((b2gd_alpha*tstop_use- &
                          dlog(cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0 &
                              + exp(-2.D0*b2gd_alpha*tstop_use) &
                              - exp(-2.D0*b2gd_alpha*tstop_use)*cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0) &
                         )/b2gd_alpha+3d0/4d0*pi) &
               + 1.D0/final_dens/R_final**5 &
                 *  exp(b2gd_alpha*tstop_use)*cos(b2gd_theta_0) &
                 / (-1.D0-sin(b2gd_theta_0) &
                    + exp(2.D0*b2gd_alpha*tstop_use)*sin(b2gd_theta_0) &
                    -exp(2.D0*b2gd_alpha*tstop_use)) &
                 * sin((-b2gd_alpha*tstop_use- &
                         dlog(cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0 &
                            + exp(-2.D0*b2gd_alpha*tstop_use) &
                            - exp(-2.D0*b2gd_alpha*tstop_use)*cos(pi/4.D0+b2gd_theta_0/2.D0)**2.D0) &
                       )/b2gd_alpha+3d0/4d0*pi) &
               )

        A_grid(itheta) = (C_sin - C_cos) / 2d0
        B_grid(itheta) = (C_sin + C_cos) / 2d0

      end do
!$OMP END DO
!$OMP END PARALLEL

      ! Compute alms.
      lsign = +1d0
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP PRIVATE(l,lsign,IA,IB) &
      !$OMP SHARED(lmax,alm,Nuse,A_grid,B_grid, handedness_sign)
      !$OMP DO SCHEDULE(static)
      do l=1, lmax
         
         ! Invert lsign.
         if (l/2==l/2.d0) then
            lsign=+1d0
         else
            lsign=-1d0
         end if

         ! Compute integrals.
         ! Use precomputed A(theta) and B(theta).
         IA = bianchi2_sky_comp_IX(l, Nuse, A_grid)
         IB = bianchi2_sky_comp_IX(l, Nuse, B_grid)

         ! Compute alm for a given 1. Only m=1 is non-zero.
         alm(l,1) = -lsign * pi * cmplx(- handedness_sign * (IB - IA), (IA + IB))

      end do
      !$OMP END DO
      !$OMP END PARALLEL

!      write(*,'(a)') ' Percent complete: 100.0%'
 
      ! Convert Delta_T/T alm computed to Delta_T alm.
      alm = BIANCHI2_CMB_T * alm

      ! Rotate alms if Euler angles present and ar least one angle non-zero.
      if(present(alpha) .and. present(beta) .and. present(gamma) &
           .and. (abs(alpha)+abs(beta)+abs(gamma) > ZERO_TOL) ) then

         write(*,*) ' Rotation of the alm '
         allocate(dl(-lmax:lmax,-lmax:lmax),stat=fail)
         allocate(alm_rotate(0:lmax,0:lmax), stat=fail)
         if(fail/=0) then
            call bianchi2_error(BIANCHI2_ERROR_MEM_ALLOC_FAIL, &
                    'bianchi_sky_rotate')
         endif
         alm_rotate = (0d0,0d0)


         do l = 1, lmax

            ! Computing the Wigner functions is the bottleneck for
            ! computation time.  If necessary, we could optimise this
            ! by precomputing dlmn(pi/2) and then using FFTs to
            ! compute dlmn(beta).
            call s2_dl_beta_operator(dl, real(beta,s2_dp), l)
  
            do m = 0,l

               ! Calculation of  Dl,m,+-1.
               Dm_p1 = exp(-icmpx*m*alpha) * dl(m, 1) * exp(-icmpx*gamma)
               Dm_m1 = exp(-icmpx*m*alpha) * dl(m,-1) * exp( icmpx*gamma)
               
               ! Rotation of the alm.
               alm_rotate(l,m) = - Dm_m1*conjg(alm(l,1)) + Dm_p1*alm(l,1)

            end do
         end do  

         ! Initialise sky object with alm_rotate.
         b%sky = s2_sky_init(alm_rotate, lmax, lmax, nside)    

         deallocate(alm_rotate)
         deallocate(dl)

      else

         write(*,'(a)') ' No rotation needed '

         ! Initialise sky object with alm.
         b%sky = s2_sky_init(alm, lmax, 1, nside)
      endif

      ! Set initialised status.
      b%init = .true.

      call bianchi2_sky_compute_map(b,nside)

      ! Free memory.
      deallocate(alm)
      deallocate(A_grid)
      deallocate(B_grid)

    end function bianchi2_sky_init_alm



    !--------------------------------------------------------------------------
    ! bianchi2_sky_free
    !
    !> Free all memory associated with a bianchi2 object.
    !!
    !!  \param[inout] b Bianchi2 object to free.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
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
    !> Write parameters of the Bianchi2 simulation to standard output.
    !!
    !! Variables:
    !!  \param[inout] b Bianchi2 object containing parameter attributes to write.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
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
    !> Compute the map of the bianchi2 sky, assuming the alms are already
    !! defined.
    !!
    !!  \param[inout] b Bianchi2 object containing alms to compute map of.
    !!  \param[in] nside Healpix nside to compute map at.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
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
    !> Compute the alm of the bianchi2 sky, assuming the map is already
    !! defined.
    !!
    !!  \param[inout] b Bianchi2 object containing sky to compute alms of.
    !!  \param[in] lmax Maximum harmonic l to consider when computing alms.
    !!  \pram[in] mmax Maximum harmonic m to consider when computing alms.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a> 
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
    !> Get sky variable from the passed bianchi2 object.
    !!
    !!   \note Initialises a new sky as a copy of the bianchi2 sky.
    !!   \note The returned sky is subsequently independed of the sky stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !!   \param[in] b Bianchi2 object containing sky to get.
    !!   \retval sky Object sky variable returned.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
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
    !> Get alms contained in a bianchi2 sky object.  Alms must already be
    !! computed.
    !!
    !!   \note Error occurs if alms are not already computed.
    !!
    !!   \param[in] b Bianchi2 object containing sky that in turn contained the alms
    !!     to get.
    !!   \retval alm Alms of simulated bianchi map extracted.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
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
    !> Apply Gaussian beam with specified FWHM.  (FWHM must be passed in
    !! arcmin.)
    !!
    !!   \note Error occurs if alms are not already computed.
    !!   
    !!   \param[inout] b Bianchi2 object containing sky that is to be comvolved with
    !!     the beam.
    !!   \param[in] fwhm Gaussian beam FWHM to use (specified in arcmin).
    !!   \param[in] lmax Maximum harmonic l to consider.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
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
    !> Write the bianchi2 simulated sky to a file. 
    !!
    !!  \param[inout] b Bianchi2 object to save sky of.
    !!  \param[in] filename Name of output file.
    !!  \param[in] file_type Type of output file, either fits map or sky file 
    !!    (see s2_sky_mod for more details).  Integer flag that may be 
    !!    either S2_SKY_FILE_TYPE_MAP or S2_SKY_FILE_TYPE_SKY.
    !!  \param[in] [comment] Optional comment to append to file header.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
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
    !> Rotate the simulated sky of the bianchi2 object.  
    !!
    !!  \param[inout] b Bianchi2 object containing sky to be rotated.
    !!  \param[in] alpha Alpha Euler angle of the rotation.
    !!  \param[in] beta Beta Euler angle of the rotation.
    !!  \param[in] gamma Gamma Euler angle of the rotation.
    !!  \param[in] lmax Maximum harmonic l to consider.
    !!  \param[in] nside Healpix nside to compute map at.
    !!  \param[in] rotation_alm Logical to specify the space used for the rotation.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
    !! \authors Thibaut Josset
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

      ! Set icmpx=sqrt(-1).
      icmpx = cmplx(0d0,1d0)

      ! Check object initialised.
      if(.not. b%init) then
        call bianchi2_error(BIANCHI2_ERROR_NOT_INIT, 'bianchi2_sky_rotate')
      end if
      

      ! Rotate if Euler angles present and at least one angle non-zero.
      if(present(alpha) .and. present(beta) .and. present(gamma) &
           .and. (abs(alpha)+abs(beta)+abs(gamma) > ZERO_TOL) ) then
         
         ! Rotation pixel by pixel or rotation of the alm.
         if (.not. rotation_alm) then

            write(*,'(a)') ' Rotation pixel by pixel with s2_sky_rotate '
            call s2_sky_rotate(b%sky, alpha, beta, gamma)

         else

            write(*,'(a)') ' Rotation of the alm '
            allocate(dl(-lmax:lmax,-lmax:lmax),stat=fail)
            allocate(alm(0:lmax,0:1), stat=fail)
            allocate(alm_rotate(0:lmax,0:lmax), stat=fail)
            if(fail/=0) then
               call bianchi2_error(BIANCHI2_ERROR_MEM_ALLOC_FAIL, &
                    'bianchi_sky_rotate')
            endif
            alm_rotate = 0e0

            ! Get alm.
            call bianchi2_sky_compute_alm(b,lmax,1)
            call bianchi2_sky_get_alm(b, alm)

            ! Perform rotation in harmonic space.
            ! Noting Bianchi alms only non zero for m=+-1).
            do l = 1,lmax

               call s2_dl_beta_operator(dl, real(beta,s2_dp), l)
            
               do m = 0,l


                  ! Calculation of  Dl,m,+-1.
                  Dm_p1 = exp(-icmpx*m*alpha) * dl(m, 1) * exp(-icmpx*gamma)
                  Dm_m1 = exp(-icmpx*m*alpha) * dl(m,-1) * exp( icmpx*gamma)

                  ! Rotation of the alm.
                  alm_rotate(l,m) = - Dm_m1*conjg(alm(l,1)) + Dm_p1*alm(l,1)
             
               end do

            end do
           
            ! Make free b%sky.
            call s2_sky_free(b%sky)

            ! Compute sky object with rotated alms.
            b%sky = s2_sky_init(alm_rotate, lmax, lmax, nside)
           
            call bianchi2_sky_compute_map(b,nside)

            ! Dellocate memory used for rotating alms.
            deallocate(alm_rotate)
            deallocate(alm)
            deallocate(dl)

         end if

       else
          write(*,'(a)') ' No rotation needed '
       endif

    end subroutine bianchi2_sky_rotate


    !--------------------------------------------------------------------------
    ! Routines to compute Bianchi induced fluctuations
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! get_results
    !
    !> This returns the values found at the end of the integral
    !!
    !!  \param[in] tstop_use
    !!  \retval R_final
    !!  \retval RH_final
    !!  \retval cos_bit_final
    !!  \retval sin_bit_final
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
    !! \authors <a href="http://www.google.co.uk/search?q=Anthony%20Lasenby">
    !!          Anthony Lasenby</a>
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
      integer :: neqn,irelab,ifail=1
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

      ! Next bit sets up initial value for R.

      vals(1)=1

      ! Next bit sets up initial value for RH.

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
    !>
    !!
    !!  \param[in] t
    !!  \param[in] vars
    !!  \retval ff
    !!
    !! \authors <a href="http://www.google.co.uk/search?q=Anthony%20Lasenby">
    !!          Anthony Lasenby</a>
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
    !>
    !!
    !!  \param[inout] b
    !!  \param[in] vals
    !! 
    !! \authors <a href="http://www.google.co.uk/search?q=Anthony%20Lasenby">
    !!          Anthony Lasenby</a>
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
    !>
    !!
    !!  \retval tau_needed
    !!
    !! \authors <a href="http://www.google.co.uk/search?q=Anthony%20Lasenby">
    !!          Anthony Lasenby</a>
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
    !>
    !!  \param[in] a
    !!
    !! \authors <a href="http://www.google.co.uk/search?q=Anthony%20Lasenby">
    !!          Anthony Lasenby</a>
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







    !--------------------------------------------------------------------------
    ! Bianchi integrals
    !--------------------------------------------------------------------------
 
    !--------------------------------------------------------------------------
    ! bianchi2_sky_comp_IX
    !
    !> Compute IA and IB integrals required in Bianchi2 simulation.
    !!
    !!   \param[in] l Harmonic l to compute IA_l or IB_l for.
    !!   \param[in] Nuse Number of terms in use in the integration.
    !!   \param[in] X_grid Contains the values A_grid(itheta) or B_grid(itheta).
    !!   \retval IX Value of the integral.
    !!
    !! \authors Thibaut Josset
    !--------------------------------------------------------------------------  

    function bianchi2_sky_comp_IX(l, Nuse, X_grid) result(IX)

      use bianchi2_plm1table_mod

      integer, intent(in) :: l, Nuse
      real(s2_dp), intent(in) :: X_grid(0:) ! X = A or B.
      real(s2_dp) :: IX
      
      real(s2_dp) :: integrand, Plm1
      integer :: itheta
      real(s2_dp) :: dtheta, theta


      IX = 0d0
      theta = 0d0
      dtheta = (pi - 0d0) / Real(Nuse, s2_dp)
      
      do itheta = 0, Nuse-1
         
         if (l>0 .and. l<=BIANCHI2_PLM1TABLE_LMAX .and. &
                Nuse == BIANCHI2_PLM1TABLE_NTHETA .and. &
                .not. BIANCHI2_DISABLE_PLM1TABLE) then

            Plm1 = bianchi2_plm1table_getval(l,theta)

         else

            Plm1 = plgndr(l, 1, cos(theta))

         end if

         integrand = sin(theta) * X_grid(itheta) * Plm1
         IX = IX + integrand*dtheta
         theta = theta + dtheta
         
      end do

      IX = IX * sqrt( (2d0*l+1d0) / real(4d0*pi*l*(l+1d0), s2_dp) ) 

      return
    end function bianchi2_sky_comp_IX
 

    !--------------------------------------------------------------------------
    ! Special function routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! plgndr
    !
    !> Computes the associated Legendre function for l and m.
    !!  Adapted from numerical recipes 
    !!
    !!   \note Numerical recipies comment:
    !!     Computes the associated Legendre polynomial P_m^l(x).
    !!
    !!   \param[in] l Legendre function l parameter.
    !!   \param[in] m Legendre function m parameter.
    !!   \param[in] x Point to evaluate specified Legendre funtion at.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
    !--------------------------------------------------------------------------

    function plgndr(l,m,x)

      integer :: l, m
      real(s2_dp) :: plgndr, x

      integer :: i, ll
      real(s2_dp) ::  fact, pll, pmm, pmmp1, somx2

      if(m<0 .or. m>l .or. abs(x)>1d0) stop 'bad arguments in plgndr'

      pmm=1.                     ! Compute Pmm.
      if(m.gt.0) then
         somx2=sqrt((1.-x)*(1.+x))
         fact=1.
         do i=1,m
            pmm=-pmm*fact*somx2
            fact=fact+2.
         enddo
      endif
      if(l.eq.m) then
         plgndr=pmm
      else
         pmmp1=x*(2*m+1)*pmm      ! Compute Pm m+1.
         if(l.eq.m+1) then
            plgndr=pmmp1
         else                    ! Compute Pm  l , l > m+ 1.
            do ll=m+2,l
               pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
               pmm=pmmp1
               pmmp1=pll
            enddo
            plgndr=pll
         endif
      endif
      return
    end function plgndr


end module bianchi2_sky_mod
