module mod_defdfpmtbn
implicit none
contains
      subroutine defdfpmtbn
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfpmtbn.
!
!     COM
!_C    Valeurs par defaut specifiees pour nbr de cycl en debut de calcul
!_C    au cours desquelles mut n'est pas mis a jour et freq en it de
!_C    mise a jour de mut.
!
!***********************************************************************
!
      use para_fige
   use kcle
   use modeleturb
implicit none
!
!-----------------------------------------------------------------------
!
      icytur0=100
      ncyturb=5
!
      kicytur0=1
      kncyturb=1
!
      return
      end
end module
