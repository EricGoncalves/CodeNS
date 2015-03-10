      subroutine def
implicit none
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut pour tout le code.
!
!***********************************************************************
!
      call defcpbd
      call defdffw
      call defdfgm
      call defdfnm
      call defdfnzst
      call defdfph
      call defdfpmcfg
      call defdfpmdtd
      call defdfpmdtg
      call defdfpmdsd
      call defdfpmimd
      call defdfpmtbn
      call defdfst
      call defdftl1
      call defintn
      call defsecpfw
!
      return
      end
