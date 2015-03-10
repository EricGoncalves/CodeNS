module mod_tcmd_cpfw
implicit none
contains
      subroutine tcmd_cpfw(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action cpfw.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
implicit none
integer :: imot
integer :: nmot
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
!
      return
      end
end module
