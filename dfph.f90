      subroutine dfph
!
!***********************************************************************
!
!     ACT
!_A    Calcul de composes de gamma (rapport des chaleurs specifiques du gaz)
!_A    pour completer la definition du gaz.
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
      use proprieteflu
implicit none
!
!-----------------------------------------------------------------------
!
      gam1=gam-1.
      gam2=.5*gam1
      gam3=1./gam
      gam4=1./gam1
      gam5=gam*gam1
!
      return
      end
