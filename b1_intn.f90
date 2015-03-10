module mod_b1_intn
implicit none
contains
      subroutine b1_intn
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'intn'.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use sortiefichier
   use kcle
   use schemanum
use mod_convich
implicit none
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
      character(len=24) ::  cnumt
!
      call convich(knumt,cnumt)
!
       form='(/,2x,''initialisation du temps "numerique"'',/' &
             //'2x,''-----------------------------------'',/' &
             //'2x,''numt                     : '',11x,i5,2x,a)'
      write(imp,form) numt,cnumt
!
      return
      end
end module
