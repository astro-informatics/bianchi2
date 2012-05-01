!------------------------------------------------------------------------------
! bianchi2_about
!
!> Display information about the BIANCHI2 package.
!!
!!   \note [-help] Display usage information.
!!
!! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!! \authors Thibaut Josset
!------------------------------------------------------------------------------

program bianchi2_about

  use s2_types_mod

  implicit none

  ! Parse input parameters.
  call parse_options()

  ! Display info.
  write(*,'(a)') "=========================================================="
  write(*,'(a)') "BIANCHI2 package to simulate bianchi-induced temperature variations"
  write(*,'(a)') "By Jason McEwen"
  write(*,'(a)') "   Anthony Lasenby"
  write(*,'(a)') "   Thibaut Josset"

  write(*,'(a)') "See www.jasonmcewen.org for more information."
  write(*,'(a)') "See LICENSE.txt for license details."

  write(*,'(a,a)') "Version: ", BIANCHI2_VERSION
  write(*,'(a,a)') "Build: ", BIANCHI2_BUILD
  write(*,'(a)') "=========================================================="


 !----------------------------------------------------------------------------

  contains


    !---------------------------------------------------------------------
    ! parse_options
    !
    !> Parse the options passed when program called.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
    !! \authors Thibaut Josset
    !---------------------------------------------------------------------

    subroutine parse_options()

      use extension, only: getArgument, nArguments
     
      implicit none
      
      integer :: n, i
      character(len=S2_STRING_LEN) :: opt
      character(len=S2_STRING_LEN) :: arg
      
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
            write(*,'(a)') 'Usage: bianchi2_about'
            stop

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program bianchi2_about




