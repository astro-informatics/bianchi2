!------------------------------------------------------------------------------
! bianchi2_about
!
!> Display information about the BIANCHI2 package.
!!
!! \section Usage Usage
!! \code 
!! bianchi2_about
!! \endcode
!!
!! \authors <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!! \authors Thibaut Josset
!------------------------------------------------------------------------------

program bianchi2_about

  use s2_types_mod

  implicit none

  ! Display info.
  write(*,'(a)') "=========================================================="
  write(*,'(a)') "BIANCHI2 package to simulate Bianchi type VIIh induced"
  write(*,'(a)') "temperature fluctuations in CMB maps when incorporating a"
  write(*,'(a)') "cosmological constant."
  write(*,'(a)') "By Jason McEwen, Anthony Lasenby and Thibaut Josset"
  write(*,'(a)') "See www.jasonmcewen.org for more information."
  write(*,'(a)') "See LICENSE.txt for license details."

  write(*,'(a,a)') "Version: ", BIANCHI2_VERSION
  write(*,'(a,a)') "Build: ", BIANCHI2_BUILD
  write(*,'(a)') "=========================================================="

end program bianchi2_about




