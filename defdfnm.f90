module mod_defdfnm
implicit none
contains
      subroutine defdfnm
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour dfnm.
!
!     COM
!_C    reservation d'un point fictif sur chaque direction dans chaque sens,
!_C    multi-grille simple (sans maillages emboites ni FMG),
!_C    calcul sur un seul niveau de grille,
!
!***********************************************************************
!
      use para_fige
      use maillage
      use kcle
      use schemanum
      use proprieteflu
implicit none
!
!-----------------------------------------------------------------------
!
      nfi=1
      knfi=1
!
      kfmg=0
      kkfmg=1
!
      lgx=1
      klgx=1
!
      kcg=1
      kkcg=1
!
      ischema=1
      kischema=1
!
      muscl=0
      kmuscl=1
      ilim=0
      kilim=1
      xk=0.
      kxk=1
!
      kprec=0
      kkprec=1
      cte=0.5
      kcte=1
      kvisq=0
      kkvisq=1
!
      ql=0.
      kql=1
      pinfl=1.
      kpinfl=1
!
      return
      end subroutine
end module
