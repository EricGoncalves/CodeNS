      subroutine b1_svgr(disc)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'svgr'.
!
!     INP
!_I    l       : arg int    ; numero de domaine
!_I    disc    : arg char   ; changement de discretisation (centre/noeud)
!_I    imp     : com int    ; unite logiq, sorties de controle
!
!***********************************************************************
!
      use sortiefichier
implicit none

!-----------------------------------------------------------------------
!
      character(len=1316) :: form
      character(len=4 ) :: disc
!
       form='(/,2x,''sauvegarde du maillage'',/' &
             //'2x,''----------------------'',/' &
             //'2x,''discretisation           : '',12x,a)'
!
      write(imp,form) disc
!
      return
      end
