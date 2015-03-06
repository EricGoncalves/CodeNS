      subroutine b1_crdms( &
                 l,ni,nj,nk)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'crdms'.
!
!     INP
!_I    td         : arg char             ; type de domaine (struct./non struct.)
!_I    ni         : arg int              ; nbr de plans indices i
!_I    nj         : arg int              ; nbr de plans indices j
!_I    nk         : arg int              ; nbr de plans indices k
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
	  use sortiefichier
!
!-----------------------------------------------------------------------
!
      character *1316 form
!
       form='(/,2x,''creation d''''un domaine '',/' &
             //'2x,''----------------------'',/' &
             //'2x,''numero                   : '',11x,i5/' &
             //'2x,''nb de pts du maillage fin: '',11x,3i5)'
      write(imp,form) l,ni,nj,nk
!
      return
      end
