module mod_b3_crdms
implicit none
contains
      subroutine b3_crdms(l,ni,nj,nk)
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
   use maillage
implicit none
integer :: l
integer :: ni
integer :: nj
integer :: nk
!
!-----------------------------------------------------------------------
!
      if(l.eq.1) then
       nptot=ni*nj*nk
      else
       nptot=nptot+ni*nj*nk
      endif
!
      if(l.eq.lzx) write(imp,999) nptot
!
 999  format(   /2x,'Nb tot de pts (hors fic.):',9x,i8, &
                /2x,'-------------------------')
!
      return
      end subroutine
end module
