      subroutine c_dfph(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfph.
!
!     INP
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    gam        : com real             ; rapport des chaleurs specifiques
!
!     OUT
!_O    gam1       : com real             ; rap chal spec -1
!_O    gam2       : com real             ; (rap chal spec -1)/2
!_O    gam3       : com real             ; 1/rap chal spec
!_O    gam4       : com real             ; 1/(rap chal spec -1)
!_O    gam5       : com real             ; rap chal spec/(rap chal spec -1)
!
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
      call tcmd_dfph(mot,imot,nmot)
!
      if(kimp.ge.1) then
        call b1_dfph
      endif
!
      call dfph
!
      return
      end
