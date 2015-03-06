      subroutine rdcmd(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Lecture d'une commande et transformation en une suite de mots.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
!
!-----------------------------------------------------------------------
!
      character *1316 command
      character *32 mot(nmx)
      dimension imot(nmx)
!
      call gtcmd(command,lgcmd)
      call splcmd(command,lgcmd,mot,imot,nmot)
!
      return
      end
