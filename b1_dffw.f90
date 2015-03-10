module mod_b1_dffw
implicit none
contains
      subroutine b1_dffw
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dffw'.
!
!     INP
!_I    equat      : com char             ; type d'equations modelisant l'ecoule-
!_I                                        ment
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    klomg      : com int              ; cle pour rotation du repere relatif
!_I    omg        : com real             ; vitesse rotation du repere relatif
!
!-----------------------------------------------------------------------
!
      use sortiefichier
   use chainecarac
   use definition
   use kcle
use mod_convich
implicit none
!
!-----------------------------------------------------------------------

      character(len=1316) :: form
      character(len=24) ::  cequat,cklomg,comg
!
      call convich(kequat,cequat)
      call convich(kklomg,cklomg)
      call convich(komg,comg)
!
       form='(/,2x,''definition de l''''ecoulement'',/' &
             //'2x,''--------------------------'',/' &
             //'2x,''equat                    : '',9x,a,2x,a/' &
             //'2x,''klomg                    : '',11x,i5,2x,a/' &
             //'2x,''omg                      : '',e12.6,''t/mn'',' &
             //'2x,a)'
      write(imp,form) equat,cequat, &
                      klomg,cklomg, &
                      omg,comg
!
      return
      end
end module
