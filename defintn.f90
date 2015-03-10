      subroutine defintn
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    intn.
!
!     COM
!_C    Par defaut le calcul debute au temps numerique 0 .
!
!***********************************************************************
!
      use para_fige
      use kcle
      use schemanum
implicit none
!
!-----------------------------------------------------------------------
!
      numt   =0
      knumt  =1
!
      return
      end
