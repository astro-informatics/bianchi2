!------------------------------------------------------------------------------
! bianchi2_lut_mod
!
!> Lookup table for associated Legendre function for case where m=1.
!! Note tailored to Bianchi2 code since input
!! for associated Legendre function is theta not x.
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
            bianchi2_lut_get_table, &
            bianchi2_lut_plgndr

  !---------------------------------------
  ! Data types
  !---------------------------------------
  type, public :: bianchi2_lut
     !> Maximal spherical harmonic to consider.
     integer :: lmax_LUT
     !> Number of points for theta.
     integer :: Nuse_LUT
     !> Contains the values of the Legendre functions.
     real(s2_dp), dimension(:,:), allocatable, public :: plgndr_table
     !> Logical to read the LUT.
     logical :: read_LUT
     !> Filename of the LUT.
     character(len=S2_STRING_LEN) :: filename_LUT
  end type bianchi2_lut

  !--------------------------------------
  ! Global Data
  !--------------------------------------
  !> Global data lut.
  type(bianchi2_lut), public :: lut

!----------------------------------------------------------------------

  contains

    !-------------------------------------------------------------------------
    ! bianchi2_lut_init
    !
    !> Initalizes the variable lut (except plgndr_table).
    !!
    !!   \param[in] read_LUT Logical to read the look_up_table.
    !!   \param[in] filename_LUT File containing the LUT.
    !!
    !! \authors Thibaut Josset
    !-------------------------------------------------------------------------
    subroutine bianchi2_lut_init(read_LUT,filename_LUT)
      
      logical :: read_LUT
      character(len=S2_STRING_LEN) :: filename_LUT

      lut%lmax_LUT = 0
      lut%Nuse_LUT = 0
      lut%read_LUT = read_LUT
      lut%filename_LUT = filename_LUT

    end subroutine bianchi2_lut_init

    !-------------------------------------------------------------------------
    ! bianchi2_lut_write
    !
    !> Writes in a file the table of Legendre functions
    !! for 1<=l<=lmax and 0<=itheta<=Nuse-1.
    !!
    !!   \param[in] lmax Maximum spherical harmonic l to consider.
    !!   \param[in] Nuse Number of terms used to compute IA and IB.
    !!                   (cf. Bianchi2_sky_mod)
    !!   \paramn[in] filename_LUT File containing the LUT.
    !!
    !! \authors Thibaut Josset
    !-------------------------------------------------------------------------
    subroutine bianchi2_lut_write(lmax,Nuse,filename_LUT)

      integer :: lmax, Nuse
      character(len=S2_STRING_LEN) :: filename_LUT

      integer :: l, itheta
      real(s2_dp) :: theta,dtheta,plm1

      open(10,file=trim(filename_LUT),form='unformatted',&
           status='replace',action='write')

      ! Write the sizes of the table.
      write(10) lmax
      write(10) Nuse

      ! Write the data.
      dtheta = (PI - 0d0) / Real(Nuse, s2_dp)
      do l=1, lmax
         theta = 0d0
         do itheta=0, Nuse-1
            ! Compute the data with plgndr.
            plm1=bianchi2_lut_plgndr(l,1,cos(theta))
            write(10) plm1
            theta = theta + dtheta
         end do
      end do

      close(10)

    end subroutine bianchi2_lut_write

    !-------------------------------------------------------------------------
    ! bianchi2_lut_get_table
    !
    !> Gets all the data in a matrix.
    !!
    !!   \param[in] lmax Maximum spherical harmonic l to consider.
    !!   \param[in] Nuse Number of terms used to compute IA and IB.
    !!
    !! \authors Thibaut Josset
    !-------------------------------------------------------------------------
    subroutine bianchi2_lut_get_table(lmax,Nuse)

      integer, intent(in) :: lmax,Nuse
      
      integer :: itheta,i,l

      ! If one wants to use table, check if it is readable.
      if (lut%read_LUT==.true.) then

         open(10,file=trim(lut%filename_LUT),form='unformatted',&
              status='old',action='read')

         read(10) lut%lmax_LUT
         read(10) lut%Nuse_LUT
         if(lmax .ne. lut%lmax_LUT) then
            call bianchi2_error(BIANCHI2_ERROR_PLM1TABLE_L_INVALID, &
                 'bianchi2_plm1table_getval')
            lut%read_LUT=.false. ! Continue without using the table.

         elseif(Nuse .ne. lut%Nuse_LUT) then
            call bianchi2_error(BIANCHI2_ERROR_PLM1TABLE_THETA_INVALID, &
                 'bianchi2_plm1table_getval')
            lut%read_LUT=.false. ! Continue without using the table.

         else
            allocate(lut%plgndr_table(1:lmax,0:Nuse-1))
            ! Transform the array in a convenient matrix.
            do l=1, lmax
               do itheta=0, Nuse-1
                  read(10) lut%plgndr_table(l,itheta)
               end do
            end do
            close(10)
         end if
      end if

    end subroutine bianchi2_lut_get_table


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
