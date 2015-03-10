module mod_c_secpfw
implicit none
contains
      subroutine c_secpfw(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action secpfw.
!
!***********************************************************************
!
      use para_fige
   use sortiefichier
use mod_tcmd_secpfw
use mod_b1_secpfw
implicit none
integer :: imot
integer :: nmot
integer :: lgr
integer :: lgrd
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
!
      dimension imot(nmx)
      dimension lgr(nobj)
!
      call tcmd_secpfw( &
                 mot,imot,nmot, &
                 lgr,lgrd)
!
      if(kimp.ge.1) then
            call b1_secpfw(lgr,lgrd)
      endif
!
      return
      end
end module
