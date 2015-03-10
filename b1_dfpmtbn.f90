      subroutine b1_dfpmtbn
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfpmtbn'.
!
!     INP
!_I    icytur0    : arg int              ; nbr de cycl en deb de calcul au cours
!_I                                        desquelles mut n'est pas mis a jour
!_I    ncyturb    : arg int              ; freq en it de mise a jour de mut
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use sortiefichier
   use modeleturb
   use kcle
implicit none
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
      character(len=24) ::  cicytur0,cncyturb
!
      call convich(kicytur0,cicytur0)
      call convich(kncyturb,cncyturb)
!
       form='(/,2x,''def. des parametres pour modele de turbulence'',/' &
             //'2x,''---------------------------------------------'',/' &
             //'2x,''icytur0                  : '',11x,i5,2x,a/' &
             //'2x,''ncyturb                  : '',11x,i5,2x,a)'
      write(imp,form) icytur0,cicytur0, &
                      ncyturb,cncyturb
!
      return
      end
