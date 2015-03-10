module mod_defdfnzst
implicit none
contains
      subroutine defdfnzst
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfnzst.
!
!***********************************************************************
!
      use kcle
      use constantes
   use definition
implicit none
!
!-----------------------------------------------------------------------
!
      ronz=reelmx
      kronz=0
!
      anz=reelmx
      kanz=0
!
      tnz=reelmx
      ktnz=0
!
      dnz=reelmx
      kdnz=0
!
      pnz=reelmx
      kpnz=0
!
      rnz=reelmx
      krnz=0
!
      return
      end
end module
