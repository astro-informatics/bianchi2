!------------------------------------------------------------------------------
! bianchi2_sim -- BIANCHI2 simulation program
!
!> Simulate a Bianchi2 VII_h model of the CMB, incorporating a cosmological
!! constant.
!!
!! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!! \authors <a href="http://www.google.co.uk/search?q=Anthony%20Lasenby">
!!          Anthony Lasenby</a>
!! \authors Thibaut Josset
!------------------------------------------------------------------------------

program bianchi2_sim

  use s2_types_mod
  use s2_sky_mod, only: S2_SKY_FILE_TYPE_MAP, S2_SKY_FILE_TYPE_SKY
  use bianchi2_sky_mod
  use bianchi2_error_mod
  use pix_tools, only: nside2npix

  use extension, only: getArgument, nArguments
  use paramfile_io, only: paramfile_handle, parse_init, parse_int, &
    parse_real, parse_double, parse_lgt, parse_string, concatnl

  implicit none

  character(len=S2_STRING_LEN) :: filename_param
  character(len=S2_STRING_LEN) :: description
  character(len=S2_STRING_LEN) :: line
  type(paramfile_handle) :: handle

  real(s2_dp), parameter :: OMEGA_MATTER_LOWER = 0d0
  real(s2_dp), parameter :: OMEGA_MATTER_UPPER = 1d0
  real(s2_dp), parameter :: OMEGA_MATTER_DEFAULT = 0.3d0
  real(s2_dp), parameter :: OMEGA_LAMBDA_LOWER = 0d0
  real(s2_dp), parameter :: OMEGA_LAMBDA_UPPER = 1d0
  real(s2_dp), parameter :: OMEGA_LAMBDA_DEFAULT = 0.6d0
  real(s2_dp), parameter :: H_LOWER = 0d0
  real(s2_dp), parameter :: H_UPPER = 10d0
  real(s2_dp), parameter :: H_DEFAULT = 1d-2
  real(s2_dp), parameter :: ZE_LOWER = 1d2
  real(s2_dp), parameter :: ZE_UPPER = 1d4
  real(s2_dp), parameter :: ZE_DEFAULT = 1d3
  real(s2_dp), parameter :: WH_DEFAULT = 1d-10
  real(s2_sp), parameter :: FWHM_DEFAULT = 330d0
  logical, parameter :: RHAND_DEFAULT = .true.
  logical, parameter :: HARMONIC_SPACE_DEFAULT = .true.
  integer, parameter :: NUSE_DEFAULT = 100
  integer, parameter :: LMAX_DEFAULT = 64
  integer, parameter :: NSIDE_DEFAULT = 128
  character(len=*), parameter :: FILE_TYPE_MAP_STR = 'map'
  character(len=*), parameter :: FILE_TYPE_SKY_STR = 'sky'

  character(len=S2_STRING_LEN) :: filename_out
  character(len=S2_STRING_LEN) :: filetype_str = FILE_TYPE_MAP_STR
  integer :: filetype = S2_SKY_FILE_TYPE_MAP

  type(bianchi2_sky) :: b
  real(s2_dp) :: omega_matter, omega_lambda, h, zE, wH
  integer :: nside, lmax
  ! Default angles in degrees for user input, but converted to radians later.
  real(s2_sp) :: alpha=0e0, beta=-90e0, gamma=0e0 
  logical :: rhand = .true.
  logical :: apply_beam = .false.
  real(s2_sp) :: fwhm = FWHM_DEFAULT

  logical :: harmonic_space
  integer :: Nuse

  logical :: rotation_alm = .false. ! Choose true only for 
                                    ! simulation in real space but 
                                    ! rotation in harmonic space.

  ! Set default parameter values.
  filename_out = 'sky.fits'
  omega_matter = OMEGA_MATTER_DEFAULT
  omega_lambda = OMEGA_LAMBDA_DEFAULT

  h = H_DEFAULT
  zE = ZE_DEFAULT
  wH = WH_DEFAULT
  nside = NSIDE_DEFAULT
  lmax = LMAX_DEFAULT
  rhand = RHAND_DEFAULT

  harmonic_space = HARMONIC_SPACE_DEFAULT
  Nuse = NUSE_DEFAULT

  write(*,'(a)') '***********************************************'
  write(*,'(a)') 'BIANCHI2 VII_h rotating universe CMB simulation'
  write(*,'(a)') 'Anthony Lasenby & Jason McEwen     October 2005'
  write(*,'(a)') 'mcewen@mrao.cam.ac.uk                          '
  write(*,'(a)') '***********************************************'


  !---------------------------------------
  ! Parse parameters
  !---------------------------------------
 
  ! Initialise file parser.
  if(nArguments() == 0) then
     filename_param = ''
  else
     if(nArguments() /= 1) then
        call bianchi2_error(BIANCHI2_ERROR_SIM_NARG, 'bianchi2_sim', &
             comment_add='Usage: bianchi2_sim [input parameter filename]')
     end if
     call getArgument(1, filename_param)
  end if
  handle = parse_init(trim(filename_param))

99 continue

  ! Get omega_matter.
  write(line,'(a,f4.1,a,f4.1,a)') '(In the range ', &
       OMEGA_MATTER_LOWER, ' <= omega_matter < ', OMEGA_MATTER_UPPER, ')'
  description = concatnl('', &
       'Enter omega_matter: ', &
       line)
1 continue
  omega_matter = parse_double(handle, 'omega_matter', &
       default=OMEGA_MATTER_DEFAULT, descr=description)
  if(omega_matter <  OMEGA_MATTER_LOWER .or. &
    omega_matter >=  OMEGA_MATTER_UPPER) then
     if(handle%interactive) goto 1
     call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
          comment_add='omega_matter invalid')
  end if

  ! Get omega_lambda.
  write(line,'(a,f4.1,a,f4.1,a)') '(In the range ', &
       OMEGA_LAMBDA_LOWER, ' <= omega_lambda < ', OMEGA_LAMBDA_UPPER, ')'
  description = concatnl('', &
       'Enter omega_lambda: ', &
       line)
2 continue
  omega_lambda = parse_double(handle, 'omega_lambda', &
       default=OMEGA_LAMBDA_DEFAULT, descr=description)
  if(omega_lambda <  OMEGA_LAMBDA_LOWER .or. &
    omega_lambda >=  OMEGA_LAMBDA_UPPER) then
     if(handle%interactive) goto 1
     call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
          comment_add='omega_lambda invalid')
  end if

  ! Check model closed.
  if((omega_matter + omega_lambda) >= 1d0) then
    if(handle%interactive) then
      write(*,*)
      write(*,'(a)') ' * Only open models are allowed *'
      goto 99
    end if
    call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
          comment_add='Only open models are allowed')
  end if

  ! Get h.
  write(line,'(a,f4.1,a,f4.1,a)') '(In the range ', &
       H_LOWER, ' <= h <= ', H_UPPER, ')'
  description = concatnl('', &
       'Enter h: ', &
       line)
3 continue
  h = parse_double(handle, 'h', &
       default=H_DEFAULT, descr=description)
  if(h <  H_LOWER .or. h >  H_UPPER) then
     if(handle%interactive) goto 3
     call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
          comment_add='h invalid')

  end if

 ! Get zE.
  write(line,'(a,e8.2,a,e8.2,a)') '(In the range ', &
       ZE_LOWER, ' <= zE <= ', ZE_UPPER, ')'
  description = concatnl('', &
       'Enter zE: ', &
       line)
4 continue
  zE = parse_double(handle, 'zE', &
       default=ZE_DEFAULT, descr=description)
  if(zE <  ZE_LOWER .or. zE >  ZE_UPPER) then
     if(handle%interactive) goto 4
     call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
          comment_add='zE invalid')

  end if

  ! Get wH.
  description = concatnl('', &
       'Enter wH: ')
5 continue
  wH = parse_double(handle, 'wH', &
       default=WH_DEFAULT, descr=description)
  if(wH < 0d0) then
     if(handle%interactive) goto 5
     call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
          comment_add='wH invalid')
  end if

  ! Get right-handedness status.
  description = concatnl('', &
       'Enter right-handedness status (logical): ')
  rhand = parse_lgt(handle, 'rhand', &
       default=RHAND_DEFAULT, descr=description)

  ! Get alpha.
  description = concatnl('', &
       'Enter alpha (degrees): ')
6 continue
  alpha = parse_real(handle, 'alpha', &
       default=alpha, descr=description)
  if(alpha < -360d0 .or. alpha > 360) then
     if(handle%interactive) goto 6
     call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
          comment_add='alpha invalid')
  end if
  alpha = alpha / 180e0 * pi

  ! Get beta.
  description = concatnl('', &
       'Enter beta (degrees): ')
7 continue
  beta = parse_real(handle, 'beta', &
       default=beta, descr=description)
  if(beta < -180d0 .or. beta > 180) then
     if(handle%interactive) goto 7
     call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
          comment_add='beta invalid')
  end if
  beta = beta / 180e0 * pi

  ! Get gamma.
  description = concatnl('', &
       'Enter gamma (degrees): ')
8 continue
  gamma = parse_real(handle, 'gamma', &
       default=gamma, descr=description)
  if(gamma < -360d0 .or. gamma > 360) then
     if(handle%interactive) goto 8
     call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
          comment_add='gamma invalid')
  end if
  gamma = gamma / 180e0 * pi

  ! Get apply_beam.
  description = concatnl('', &
    'Enter apply_beam status (logical): ')
  apply_beam = parse_lgt(handle, 'apply_beam', &
    default=apply_beam, descr=description)

  if(apply_beam) then
    ! Get beam fwhm.
    description = concatnl('', &
       'Enter beam fwhm (arcmin): ')
9  continue
    fwhm = parse_real(handle, 'fwhm', &
        default=FWHM_DEFAULT, descr=description)
    if(fwhm <= 0d0) then
      if(handle%interactive) goto 9
      call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
            comment_add='fwhm invalid')
    end if
  end if

  ! Get lmax.
     description = concatnl("", &
          "Enter the maximum harmonic l (lmax) for the simulated sky: ")
10   continue
     lmax = parse_int(handle, 'lmax', &
          default=LMAX_DEFAULT, descr=description)
     if(lmax < 0) then
        if(handle%interactive) goto 10
        call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
             comment_add='lmax invalid')
     endif


  ! Get nside.
  description = concatnl("", &
       "Enter the resolution parameter (nside) for the simulated sky: ", &
       "(npix = 12*nside**2, where nside must be a power of 2)")
11 continue
  nside = parse_int(handle, 'nside', &
       default=NSIDE_DEFAULT, descr=description)
  if(nside2npix(nside) < 0) then
     if(handle%interactive) goto 11
     call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
          comment_add='nside invalid')
  endif
  
  ! Get harmonic_space.
  description = concatnl('', &
       'Enter harmonic_space status (logical): ')
  harmonic_space = parse_lgt(handle, 'harmonic_space', &
       default=HARMONIC_SPACE_DEFAULT, descr=description)


  ! Get Nuse.
  description = concatnl("", &
       "Enter the number of terms (Nuse) used for the integrations: ")
  Nuse = parse_int(handle, 'Nuse', &
       default=NUSE_DEFAULT, descr=description)
 

  ! Get filename_out.
  description = concatnl('', &
       'Enter filename_out: ')
  filename_out = parse_string(handle, 'filename_out', &
       default=trim(filename_out), descr=description)

  ! Get output file type: map or sky.
  description = concatnl('', &
       'Enter output file type (filetype={map; sky}): ')
12  continue
  filetype_str = parse_string(handle, 'filetype', &
       default=trim(filetype_str), descr=description)

  ! Set filetype integer status.
  select case (filetype_str)
     
    case (FILE_TYPE_MAP_STR)
      filetype = S2_SKY_FILE_TYPE_MAP

    case (FILE_TYPE_SKY_STR)
      filetype = S2_SKY_FILE_TYPE_SKY

    case default
       if(handle%interactive) goto 12
       call bianchi2_error(BIANCHI2_ERROR_SIM_PARAM_INVALID, 'bianchi2_sim', &
         comment_add='Invalid output file type')

  end select


  !---------------------------------------
  ! Run simulation and save sky
  !---------------------------------------

  ! Simulated bianchi2 sky.
  write(*,'(a)')


  ! Initialise bianchi2 object.
  if (harmonic_space == .true.) then
     
     write(*,'(a)') 'Computing BIANCHI2 simulation in harmonic space...'
     b = bianchi2_sky_init_alm(omega_matter, omega_lambda, h, zE, wH, rhand, &
          nside, lmax, Nuse, alpha, beta, gamma)
     write(*,'(a)') 'Simulation complete'
     write(*,'(a)')
  
  else

     write(*,'(a)') 'Computing BIANCHI2 simulation in real space...'
     b = bianchi2_sky_init(omega_matter, omega_lambda, h, zE, wH, rhand, &
          nside)

     write(*,'(a)') 'Simulation complete'
     write(*,'(a)')
  
     ! Perform rotation
     call bianchi2_sky_rotate(b, alpha, beta, gamma, lmax, nside,rotation_alm)

  end if

  ! Apply beam if required.
  if(apply_beam) then
      
    ! Recompute alms if necessary (if not necessary does nothing).
    ! Note, if compute alms directly but then perform
    ! bianchi2_sky_rotate, alms are removed and must be recomputed.
    call bianchi2_sky_compute_alm(b, lmax, lmax)

    ! If map already present then apply beam will recompute the map.
    call bianchi2_sky_apply_beam(b, fwhm, lmax)

  end if
 
  ! Save sky.
  call bianchi2_sky_write(b, filename_out, filetype)
  write(*,'(a,a)') 'Simulated map written to ', trim(filename_out)
  write(*,'(a)')

  ! Write bianchi2 simulation variables to standard output.
  write(*,'(a)') 'Simulation parameters:'
  call bianchi2_sky_param_write(b)

  ! Free bianchi2 object.
  call bianchi2_sky_free(b)

end program bianchi2_sim
