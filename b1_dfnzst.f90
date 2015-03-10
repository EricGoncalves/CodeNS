      subroutine b1_dfnzst(nonzst)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfnzst'.
!
!     INP
!_I    imp    : com int   ; unite logiq, sorties de controle
!_I    tnz    : com real  ; etat pour adimensionnement temperature
!_I    ronz   : com real  ; etat pour adimensionnement,masse volumique
!_I    anz    : com real  ; etat pour adimensionnement vitesse du son d'arret
!_I    dnz    : com real  ; etat pour adimensionnement longueur
!
!***********************************************************************
!
      use sortiefichier
implicit none
integer :: nonzst
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
!
      form='(//,2x,''definition de la normalisation''/' &
             //'2x,''------------------------------'',/' &
             //'2x,''numero de l''''etat         : '',9x,i2)'
      write(imp,form)nonzst
!
      return
      end
