      subroutine c_dftl1(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dftl1.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use sortiefichier
!
!-----------------------------------------------------------------------
!
      character *32 mot(nmx)
      dimension imot(nmx)
!
      call tcmd_dftl1(mot,imot,nmot)
!
!      if (kimp.ge.1) then
!            call b1_dftl1
!      endif
!
      return
      end
