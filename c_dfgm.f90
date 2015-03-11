module mod_c_dfgm
implicit none
contains
      subroutine c_dfgm(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfgm.
!
!     INP
!_I    config     : com char             ; type de config geometrique du calcul
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!
!     OUT
!_O    ptrans     : com real             ; distance pour periodicite
!_O    protat     : com real             ; angle(rad) pour periodicite
!
!     I/O
!_/    perio      : com real             ; periodicite geometrique en angle ou
!_/                                        distance selon config
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use sortiefichier 
use mod_tcmd_dfgm
use mod_b1_dfgm

use mod_dfgm

implicit none
integer :: imot
integer :: nmot
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
!
      call tcmd_dfgm(mot,imot,nmot)
!
      if (kimp.ge.1) then
        call b1_dfgm
      endif
!
      call dfgm
!
      return
      end subroutine
end module
