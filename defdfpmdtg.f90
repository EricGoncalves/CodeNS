      subroutine defdfpmdtg
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfpmdtg.
!
!     COM
!_C    Pas de valeurs par defaut pour la cle d'utilisation du pas de
!_C    temps local, la frequence de calcul du pas de temps, le pas de
!_C    temps constant eventuellement impose.
!_C    Par defaut le pas de temps est calcule a chaque iteration pendant
!_C    10 iterations avant d'etre calcule selon la frequence specifiee.
!
!***********************************************************************
!
      use para_fige
      use kcle
      use constantes
   use schemanum
!
!-----------------------------------------------------------------------
!
      kdtl   =intmx
      icychr0=10
      ncychro=intmx
      dt1min =reelmx
!
      kkdtl   =0
      kicychr0=1
      kncychro=0
      kdt1min =0
!
      return
      end
