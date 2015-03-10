module mod_b2_crbds
implicit none
contains
      subroutine b2_crbds(mfbe)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour les commande 'crbds'.
!
!     INP
!_I    mmb        : arg int    ; nombre de pts d'une frontiere
!_I    mpb        : arg int    ; pointeur fin de front precedente
!_I                              dans tableaux de base des front.
!_I    imp        : com int     ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use sortiefichier
   use maillage
   use boundary
implicit none
integer :: mfbe
integer :: img
integer :: mfbi
integer :: mfbim
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
!
      mfbi=nfei(mfbe)
!
      form='(/ 2x,''numero de la grille      : '',11x,i5/' &
            //'2x,''nb de pts de la frontiere: '',11x,i5/' &
            //'2x,''tot. pts des front preced: '',11x,i5)'
!
      do img=1,lgx
      mfbim=mfbi+(img-1)*mtb
      write(imp,form) img,mmb(mfbim),mpb(mfbim)
      enddo
!
      return
      end
end module
