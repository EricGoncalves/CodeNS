      subroutine c_end(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action end.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
implicit none
integer :: imot
integer :: nmot
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
!
      call b1_end
      call eend
!
      return
      end
