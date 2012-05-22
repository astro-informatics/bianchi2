!----------------------------------------------------------------------------
! bianchi2_lut_gen -- BIANCHI2 writting LUT program
!
!> Generate a Look-Up-Table containing the values of the Legendre functions
!! for a given couple (lmax,Nuse).
!!
!! \authors Thibaut Josset
!----------------------------------------------------------------------------

program bianchi2_lut_gen

  use s2_types_mod
  use bianchi2_lut_mod, only: bianchi2_lut_write
  use extension, only: getArgument, nArguments

  implicit none

  integer :: lmax, Nuse
  character(len=S2_STRING_LEN) :: filename_LUT

  integer :: n, i
  character(len=S2_STRING_LEN) :: arg
  character(len=S2_STRING_LEN) :: opt

  !------------------------------------------------------------------------
  ! Parse the arguments.
  !------------------------------------------------------------------------
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
        write(*,'(a)') '                        [-Nuse Nuse]'
        write(*,'(a)') '                        [-out filename_LUT(without extension)]'
        stop
          
     case ('-lmax')
        read(arg,*) lmax

     case ('-Nuse')  
        read(arg,*) Nuse

     case ('-out')
        filename_LUT = trim(arg)

     case default
        print '("Unknown option ",a," ignored")', trim(opt)            

     end select
  end do

  !-------------------------------------------------------------------------
  ! Write the LUT.
  !-------------------------------------------------------------------------
  call bianchi2_lut_write(lmax, Nuse, filename_LUT)
  write(*,'(a)') 'LUT written in '//trim(filename_LUT)//'.dat'

end program bianchi2_lut_gen
