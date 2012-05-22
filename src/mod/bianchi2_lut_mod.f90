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
  use paramfile_io, only: paramfile_handle, parse_init, parse_int, &
    parse_real, parse_double, parse_lgt, parse_string, concatnl

  implicit none

  private

  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: bianchi2_lut_write, &
            bianchi2_lut_check_sizes, &
            bianchi2_lut_get_table, &
            bianchi2_lut_plgndr


  contains

    !-------------------------------------------------------------------------
    ! bianchi2_lut_write
    !
    !> Writes in a file the table of Legendre functions
    !! for 1<=l<=lmax and 0<=itheta<=Nuse-1.
    !!
    !!   \param[in] lmax Maximum spherical harmonic l to consider.
    !!   \param[in] Nuse Number of terms used to compute IA and IB.
    !!                   (cf. Bianchi2_sky_mod)
    !!   \paramn[in] filename_LUT File to write the LUT.
    !!
    !! \authors Thibaut Josset
    !-------------------------------------------------------------------------
    subroutine bianchi2_lut_write(lmax,Nuse,filename_LUT)

      integer :: lmax, Nuse
      character(len=S2_STRING_LEN) :: filename_LUT

      integer :: l, itheta
      real(s2_dp) :: theta, plm1

      ! Write the sizes of the table in an independant file (easier to read).
      open(10,file=trim(filename_LUT)//'_sizes.dat')
      write(10,*) 'BIANCHI2_PLM1TABLE_LMAX =', LMAX
      write(10,*) 'BIANCHI2_PLM1TABLE_NUSE =', Nuse
      close(10)
     
      ! Write the data.
      open(10,file=trim(filename_LUT)//'.dat')
      do l=1, lmax
         do itheta=0, Nuse-1
            theta=itheta*PI/Nuse
            ! Comput the data with plgndr.
            plm1=bianchi2_lut_plgndr(l,1,cos(theta))
            write(10,*) plm1
         end do
      end do
      close(10)

    end subroutine bianchi2_lut_write


    !----------------------------------------------------------------------
    ! bianchi2_lut_check_sizes
    !
    !> Checks the sizes of the table before read it.
    !!
    !!   \param[in] lmax Maximum spherical harmonic l to consider.
    !!   \param[in] Nuse Number of terms used to compute IA and IB.
    !!   \param[inout] read_LUT_plgndr Logical to say whether
    !!                 the LUT should be used and whether it is readable.
    !!   \param[in] filename_LUT File to read the LUT.
    !!
    !! \authors Thibaut Josset
    !----------------------------------------------------------------------
    subroutine bianchi2_lut_check_sizes(lmax,Nuse,read_LUT,filename_LUT)
      
      integer, intent(in) :: lmax, Nuse
      logical, intent(inout) :: read_LUT
      character(len=S2_STRING_LEN) :: filename_LUT

      type(paramfile_handle) :: handle

      integer :: BIANCHI2_PLM1TABLE_LMAX, BIANCHI2_PLM1TABLE_NUSE

      if (read_LUT==.true.) then

         handle = parse_init(trim(filename_LUT)//'_sizes.dat')

         ! Get BIANCHI2_PLM1TABLE_LMAX
         BIANCHI2_PLM1TABLE_LMAX = parse_int(handle, 'BIANCHI2_PLM1TABLE_LMAX', &
              default=-1)
         ! Check lmax in a valid range.
         if(lmax .ne. BIANCHI2_PLM1TABLE_LMAX) then
            call bianchi2_error(BIANCHI2_ERROR_PLM1TABLE_L_INVALID, &
                 'bianchi2_plm1table_getval')
            read_LUT=.false. ! Continue without using the table.
         end if

         ! Get BIANCHI2_PLM1TABLE_NUSE
         BIANCHI2_PLM1TABLE_NUSE = parse_int(handle, 'BIANCHI2_PLM1TABLE_NUSE', &
              default=-1)
         ! Check Nuse in a valid range.
         if(Nuse .ne. BIANCHI2_PLM1TABLE_NUSE) then
            call bianchi2_error(BIANCHI2_ERROR_PLM1TABLE_THETA_INVALID, &
                 'bianchi2_plm1table_getval')
            read_LUT=.false. ! Continue without using the table.
         end if

      end if

    end subroutine bianchi2_lut_check_sizes


    !-----------------------------------------------------------------------------
    ! bianchi2_lut_get_table
    !
    !> Gets all the data in a matrix.
    !!
    !!   \param[in] lmax Maximum spherical harmonic l to consider.
    !!   \param[in] Nuse Number of terms used to compute IA and IB.
    !!   \param[inout] plgndr_table Contains the values of the Legendre functions.
    !!   \param[in] filename_LUT File to read the LUT.
    !!
    !! \authors Thibaut Josset
    !-----------------------------------------------------------------------------
    subroutine bianchi2_lut_get_table(plgndr_table,lmax,Nuse,filename_LUT)

      integer, intent(in) :: lmax, Nuse
      real(s2_dp), dimension(1:lmax, 0:Nuse-1), intent(inout) :: plgndr_table
      character(len=S2_STRING_LEN) :: filename_LUT

      real(s2_dp), allocatable :: data_table(:)
      integer :: l,itheta,i

      open(10,file=trim(filename_LUT)//'.dat')

      ! Put all the data in an array.
      allocate(data_table(0:lmax*Nuse-1))
      read(10,*)(data_table(i),i=0,lmax*Nuse-1)

      ! Transform the array in a convenient matrix.
      do l=1, lmax
         do itheta=0, Nuse-1
            plgndr_table(l,itheta)=data_table((l-1)*Nuse+itheta)
         end do
      end do

      close(10)
      deallocate(data_table)

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
