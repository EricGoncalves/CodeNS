      subroutine c_dfst(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Defintion d'un etat aerodynamique par la donnees de conditions
!_A    generatrices (Pi et Ti) et d'un nombre de Mach
!
!     INP
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!
!     OUT
!_O    pist       : com real             ; pression generatrice
!_O    tist       : com real             ; temperature generatrice
!_O    mast       : com real             ; nombre de Mach
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
      call tcmd_dfst( &
                 mot,imot,nmot, &
                 nst)
!
      if (kimp.ge.1) then
            call b1_dfst(nst)
      endif
!
      return
      end
