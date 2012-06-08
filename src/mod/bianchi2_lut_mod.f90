!------------------------------------------------------------------------------
! bianchi2_lut_mod
!
!> Lookup table for associated Legendre function for case where m=1.
!! Also contains functionality to compute associated Legendre
!! functions if a LUT is not used.
!!
!! \note Tailored to Bianchi2 code since input for associated Legendre
!! function is theta not x.
!!
!! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!! \authors Thibaut Josset
!------------------------------------------------------------------------------

module bianchi2_lut_mod

  use s2_types_mod
  use bianchi2_error_mod

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: bianchi2_lut_init, &
            bianchi2_lut_write, &
            bianchi2_lut_read, &
            bianchi2_lut_access, &
            bianchi2_lut_free, &
            bianchi2_lut_plgndr


  !---------------------------------------
  ! Data types
  !---------------------------------------

  !> Look-up-table that contains precomputed associated Legendre functions.
  type, public :: bianchi2_lut
     private
     !> Maximal spherical harmonic to consider.
     integer :: lmax_LUT = 0
     !> Number of points for theta.
     integer :: Ntheta_LUT = 0
     !> Contains the values of the Legendre functions.  
     !! Indexed by table_LUT(1:lmax,0:Ntheta-1).
     real(s2_dp), dimension(:,:), allocatable :: table_LUT
  end type bianchi2_lut


  !---------------------------------------------------------------------------

  contains


    !-------------------------------------------------------------------------
    ! bianchi2_lut_init
    !
    !> Initalizes a lut using lmax and Ntheta.
    !!
    !!   \param[in] lmax Maximal harmonic to consider.
    !!   \param[in] Ntheta Number of points for theta.
    !!   \retval lut Initialized bianchi2_lut object.
    !!
    !! \authors Thibaut Josset
    !-------------------------------------------------------------------------

    function bianchi2_lut_init(lmax,Ntheta) result(lut)
      
      integer, intent(in) :: lmax, Ntheta
      type(bianchi2_lut) :: lut

      integer :: l, itheta
      real(s2_dp) :: theta,dtheta

      lut%lmax_LUT = lmax
      lut%Ntheta_LUT = Ntheta

      allocate(lut%table_LUT(1:lmax,0:Ntheta-1))

      dtheta = (PI - 0d0) / Real(Ntheta, s2_dp)
      do l=1, lmax
         do itheta=0, Ntheta-1
            theta = itheta*dtheta
            ! Compute the data with plgndr.
            lut%table_LUT(l,itheta) = bianchi2_lut_plgndr(l,1,cos(theta))
            theta = theta + dtheta
         end do
      end do

    end function bianchi2_lut_init


    !-------------------------------------------------------------------------
    ! bianchi2_lut_write
    !
    !> Writes a lut in a binary file.
    !!
    !!   \param[in] lut Initialized bianchi2_lut object
    !!   \paramn[in] filename_LUT Name of binary file containing the LUT.
    !!
    !! \authors Thibaut Josset
    !-------------------------------------------------------------------------

    subroutine bianchi2_lut_write(lut,filename_LUT)

      type(bianchi2_lut), intent(in) :: lut
      character(len=S2_STRING_LEN), intent(in) :: filename_LUT

      integer :: l,itheta

      open(10,file=trim(filename_LUT),form='unformatted',&
           status='replace',action='write')

      ! Write the sizes of the table.
      write(10) lut%lmax_LUT
      write(10) lut%Ntheta_LUT

      ! Write the table
      do l=1, lut%lmax_LUT
         do itheta=0, (lut%Ntheta_LUT - 1)
            write(10) lut%table_LUT(l,itheta)
         end do
      end do

      close(10)

    end subroutine bianchi2_lut_write


    !-------------------------------------------------------------------------
    ! bianchi2_lut_read
    !
    !> Reads a lut from a file.
    !!
    !!   \param[in] lmax Maximum spherical harmonic l to consider.
    !!   \param[in] Ntheta Number of terms used to compute IA and IB.
    !!   \param[in] filename_LUT Name of binary file containing the LUT.
    !!   \retval lut bianchi2_lut object to get.
    !!
    !! \authors Thibaut Josset
    !-------------------------------------------------------------------------

    function bianchi2_lut_read(lmax,Ntheta,filename_LUT) result(lut)
      
      integer, intent(in) :: lmax,Ntheta
      character(len=S2_STRING_LEN) :: filename_LUT
      type(bianchi2_lut) :: lut

      integer :: l,itheta

      open(10,file=trim(filename_LUT),form='unformatted',&
              status='old',action='read')

      read(10) lut%lmax_LUT
      read(10) lut%Ntheta_LUT
      allocate(lut%table_LUT(1:lmax,0:Ntheta-1))

      if(lmax .ne. lut%lmax_LUT) then
         call bianchi2_error(BIANCHI2_ERROR_PLM1TABLE_L_INVALID, &
                 'bianchi2_plm1table_getval')
  
      elseif(Ntheta .ne. lut%Ntheta_LUT) then
         call bianchi2_error(BIANCHI2_ERROR_PLM1TABLE_THETA_INVALID, &
              'bianchi2_plm1table_getval')

      else     
         do l=1, lmax      
            do itheta=0, Ntheta-1
               read(10) lut%table_LUT(l,itheta)
            end do
         end do
      end if

      close(10)
      
    end function bianchi2_lut_read


    !--------------------------------------------------------------------------
    ! bianchi2_lut_access
    !
    !> Gives the value in the lut corresponding to a given (l,itheta).
    !!
    !!   \param[in] lut Initialized bianchi2_lut object.
    !!   \param[in] l
    !!   \param[in] itheta 
    !!   \retval plm1 Value of the Lengendre functions contained in the lut.
    !!
    !! \authors Thibaut Josset
    !--------------------------------------------------------------------------

    function bianchi2_lut_access(lut,l,itheta) result(plm1)
      
      type(bianchi2_lut), intent(in) :: lut
      integer, intent(in) :: l, itheta
      real(s2_dp) :: plm1

      plm1 = lut%table_LUT(l,itheta)

    end function bianchi2_lut_access

    !--------------------------------------------------------------------------
    ! bianchi2_lut_free
    !
    !> Deallocate the Look_Up_Table.
    !!
    !!   \param[in] lut Bianchi2_lut object to make free.
    !!
    !! \authors Thibaut Josset
    !--------------------------------------------------------------------------

    subroutine bianchi2_lut_free(lut)
      
      type(bianchi2_lut), intent(inout) :: lut

      if(allocated(lut%table_LUT)) deallocate(lut%table_LUT)
      lut%lmax_LUT = 0
      lut%Ntheta_LUT = 0

    end subroutine bianchi2_lut_free


    !--------------------------------------------------------------------------
    ! bianchi2_lut_plgndr
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

    function bianchi2_lut_plgndr(l,m,x) result(plgndr)

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

    end function bianchi2_lut_plgndr


end module bianchi2_lut_mod
