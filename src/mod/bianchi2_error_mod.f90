!------------------------------------------------------------------------------
! bianchi2_error_mod  -- BIANCHI2 library error class
!
!> Functionality to handle errors that may occur in the bianchi2 library. 
!! Public bianchi2 error codes are defined, with corresponding private error 
!! comments and default halt execution status.
!!
!! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!------------------------------------------------------------------------------

module bianchi2_error_mod

  use s2_types_mod, only: S2_STRING_LEN

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: bianchi2_error


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  integer, parameter :: BIANCHI2_ERROR_NUM = 9

  integer, public, parameter :: &
    BIANCHI2_ERROR_NONE = 0, &
    BIANCHI2_ERROR_INIT = 1, &
    BIANCHI2_ERROR_NOT_INIT = 2, &
    BIANCHI2_ERROR_INIT_FAIL = 3, &
    BIANCHI2_ERROR_MEM_ALLOC_FAIL = 4, &
    BIANCHI2_ERROR_SIM_PARAM_INVALID = 5, &
    BIANCHI2_ERROR_SIM_NARG = 6, &
    BIANCHI2_ERROR_SKY_NUM_FAIL = 7, &
    BIANCHI2_ERROR_TMPLFIT_FAIL = 8

  !> Each element of the error_comment array must have the same length, thus
  ! space with trailing space characters.  When come to use trim to remove 
  ! trailing spaces.
  ! Comment associated with each error type.
  character(len=S2_STRING_LEN), parameter :: &
    error_comment(BIANCHI2_ERROR_NUM) = &
      (/ & 
      'No error                                                                 ', &
      'Attempt to initialise object that has already been initialised           ', &
      'Object not initialised                                                   ', &
      'Object initialisation failed                                             ', &
      'Memory allocation failed                                                 ', &
      'Invalid simulation parameter                                             ', &
      'Invalid number of command line parameters                                ', &
      'Numerical routine failed                                                 ', &
      'Template fitting failed                                                  ' &
      /) 
  
  !> Default program halt status of each error type.
  logical, parameter :: &
    halt_default(BIANCHI2_ERROR_NUM) = &
      (/ &
      .false., &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  & 
      .true.,  &
      .true.,  &
      .true.,  &
      .true.  /)
  
  
  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! bianchi2_error
    !
    !> Display error message corresponding to error_code and halt program 
    !! execution if required.
    !!
    !! Variables:
    !!   \param[in] error_code Integer error code.
    !!   \apram[in] procedure Procedure name where bianchi2_error called from.  Displayed 
    !!     when error message printed to screen.
    !!   \param[in] comment_add If present, additional comment to append to default 
    !!     error comment.
    !!   \param[inout] comment_out If present the error comment is copied to comment_out
    !!     on output.
    !!   \param[in] halt_in  If present overrides default halt value.
    !!
    !! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
    !--------------------------------------------------------------------------

    subroutine bianchi2_error(error_code, procedure, comment_add, &
      comment_out, halt_in)

      integer, intent(in) :: error_code
      character(len=*), intent(in), optional :: procedure, comment_add
      character(len=*), intent(inout), optional :: comment_out
      logical, intent(in), optional :: halt_in

      logical :: halt
      character(len=*), parameter :: comment_prefix = 'BIANCHI2_ERROR: '

      !---------------------------------------
      ! Display error message
      !---------------------------------------

      if(present(procedure)) then

        if(present(comment_add)) then
	  write(*,'(a,a,a,a,a,a,a,a)') comment_prefix, 'Error ''', &
            trim(error_comment(error_code+1)), &
            ''' occured in procedure ''', &
            trim(procedure), &
            '''', &
            ' - ', trim(comment_add)
        else
          write(*,'(a,a,a,a,a,a)') comment_prefix, 'Error ''', &
            trim(error_comment(error_code+1)), &
            ''' occured in procedure ''', &
            trim(procedure), &
            ''''
        end if
 
     else

        if(present(comment_add)) then
          write(*,'(a,a,a,a)') comment_prefix, &
            trim(error_comment(error_code+1)), &
            ' - ', trim(comment_add)
        else
          write(*,'(a,a)') comment_prefix, trim(error_comment(error_code+1))
        end if

      end if

      !> Copy error comment if comment_out present.
      if(present(comment_out)) comment_out = error_comment(error_code+1)

      !---------------------------------------
      ! Halt program execution if required
      !---------------------------------------
      
      if( present(halt_in) ) then
        halt = halt_in
      else
        halt = halt_default(error_code+1)
      end if

      if( halt ) then
        write(*,'(a,a,a,a,a)') comment_prefix, &
          '  Halting program execution ', &
          'due to error ''', trim(error_comment(error_code+1)), ''''
        stop
      end if

    end subroutine bianchi2_error


end module bianchi2_error_mod
