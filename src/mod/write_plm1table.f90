program write_plm1table

implicit none


integer, parameter :: LMAX= 1024, NTHETA=1000
integer :: l, itheta
real(SELECTED_REAL_KIND(12,200)) :: theta, plm1
real(SELECTED_REAL_KIND(12,200)) :: plgndr

!! PI definition.
  real(SELECTED_REAL_KIND(12,200)), parameter :: PI = 3.141592653589793238462643383279502884197


open(1,file='newtable.f90')

do l=1, LMAX

   do itheta=0, NTHETA-1
      
      theta=itheta*PI/NTHETA
      plm1=plgndr(l,1,cos(theta))
      write(1,*) 'data PLM1_TABLE(',l,',',itheta,') / ',plm1,'/'

   end do

end do

close(1)

end program write_plm1table




function plgndr(l,m,x)

implicit none

      integer :: l, m
      real(SELECTED_REAL_KIND(12,200)) :: plgndr, x

      integer :: i, ll
      real(SELECTED_REAL_KIND(12,200)) ::  fact, pll, pmm, pmmp1, somx2

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
