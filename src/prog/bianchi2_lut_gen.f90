!----------------------------------------------------------------------------
! bianchi2_lut_gen -- BIANCHI2 writting LUT program
!
!> Generate a Look-Up-Table containing the values of the Legendre functions
!! for a given couple (lmax,Ntheta).
!!
!! \section Usage Usage
!! \code 
!! bianchi2_lut_gen [-lmax lmax]
!!                  [-Ntheta Ntheta]
!!                  [-out filename_LUT]
!! \endcode
!!
!! \param[in] lmax Maximum spherical harmonic order to consider.
!! \param[in] Ntheta Number of theta grid points to consider.
!! \param[in] out Filename of the output LUT.
!!
!! \authors Thibaut Josset
!----------------------------------------------------------------------------

program bianchi2_lut_gen

  use s2_types_mod
  use bianchi2_lut_mod

  implicit none

  type(bianchi2_lut) :: lut
  integer :: lmax, Ntheta
  character(len=S2_STRING_LEN) :: filename_LUT

  !Get the arguments.
  call parse_options

  ! Initialize the lut.
  lut = bianchi2_lut_init(lmax,Ntheta)

  ! Write the lut into a file.
  call bianchi2_lut_write(lut,filename_LUT)
  write(*,'(a)') 'LUT written in '//trim(filename_LUT)


  contains

  !------------------------------------------------------------------------
  ! Parse options
  !
  !> Parse the options passed when program called.
  !!
  !! \authors Thibaut Josset
  !------------------------------------------------------------------------

    subroutine parse_options

      use extension, only: getArgument, nArguments
     
      implicit none

      integer :: n, i
      character(len=S2_STRING_LEN) :: arg
      character(len=S2_STRING_LEN) :: opt
      
      n = nArguments()
     
      do i = 1,n,2
        
         call getArgument(i,opt)
     
         if (i == n .and. trim(opt) /= '-help') then
            write(*,'(a,a,a)') 'Option ', trim(opt), ' has no argument'
            stop
         end if
     
         if(trim(opt) /= '-help') call getArgument(i+1,arg)
         
         ! Read each argument in turn
         select case (trim(opt))
            
         case ('-help')
            write(*,'(a)') 'Usage: bianchi2_lut_gen [-lmax lmax]'
            write(*,'(a)') '                        [-Ntheta Ntheta]'
            write(*,'(a)') '                        [-out filename_LUT]'
            stop
            
         case ('-lmax')
            read(arg,*) lmax

         case ('-Ntheta')  
            read(arg,*) Ntheta
            
         case ('-out')
            filename_LUT = trim(arg)

         case default
            print '("Unknown option ",a," ignored")', trim(opt)            
            
         end select
      end do

    end subroutine parse_options
    

end program bianchi2_lut_gen
