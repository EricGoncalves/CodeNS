module mod_defsecpfw
implicit none
contains
      subroutine defsecpfw
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    secpfw.
!
!     COM
!_C    Par defaut verification de la metrique en debut de calcul,
!_C    frequence de sortie des residus specifiee,
!_C    frequence de sauvegarde et d'exploitation des resultats en cours
!_C    de calcul quasiment nulle,
!_C    sauvegarde des resultats aux centres des cellules.
!_C    Pas de valeur par defaut pour le nombre de cycles a effectuer.
!
!***********************************************************************
!
      use para_fige
      use constantes
      use kcle
      use schemanum
      use maillage
implicit none
integer :: ng
!
!-----------------------------------------------------------------------
!
      kvn   =0
      kkvn  =1
!
      ncyresi =25
      ncysave =1000000
      ncyexpl =1000000
      kncyresi=1
      kncysave=1
      kncyexpl=1
!
      discsv ='cccc'
      kdiscsv=1
!
      do ng=1,lg
      ncycle(ng)=intmx
      kncycle(ng)=0
      enddo
!
      return
      end
end module
