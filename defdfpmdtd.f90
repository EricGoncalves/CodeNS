      subroutine defdfpmdtd
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfpmdtd.
!
!***********************************************************************
!
      use para_fige
   use schemanum
   use kcle
implicit none
integer :: l
!
!-----------------------------------------------------------------------
!
      do l=1,lz
      eta(l) =1.0
      enddo
!
      do l=1,lz
      keta(l)=1
      enddo
!
      return
      end
