module mod_defdfpmdsd
implicit none
contains
      subroutine defdfpmdsd
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfpmdsd.
!
!     COM
!_C    Des valeurs par defaut sont donnees a tous les parametres
!_C    d'utilisation de la dissipation artificielle.
!
!***********************************************************************
!
      use para_fige
   use kcle
   use schemanum
implicit none
integer :: l
!
!-----------------------------------------------------------------------
!
      do l=1,lz
       ki2(l)=0.0
       ki4(l)=0.032
      enddo
!
      do l=1,lz
       kki2(l)=1
       kki4(l)=1
      enddo
!
      return
      end
end module
