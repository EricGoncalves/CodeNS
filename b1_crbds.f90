      subroutine b1_crbds( &
                 mfb,kini,l,imin,imax,jmin,jmax,kmin,kmax, &
                 indfl)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'crbds'.
!
!     INP
!_I    mfb        : arg int              ; numero de frontiere
!_I    td         : arg char             ; type de domaine (struct./non struct.)
!_I    l          : arg int              ; numero de domaine
!_I    imin       : arg int              ; indice min en i
!_I    imax       : arg int              ; indice max en i
!_I    jmin       : arg int              ; indice min en j
!_I    jmax       : arg int              ; indice max en j
!_I    kmin       : arg int              ; indice min en k
!_I    kmax       : arg int              ; indice max en k
!_I    indfl      : arg char             ; type de plan de la frontiere
!_I    kini       : arg int              ; cle creation de frontiere
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use sortiefichier
implicit none
integer :: mfb
integer :: kini
integer :: l
integer :: imin
integer :: imax
integer :: jmin
integer :: jmax
integer :: kmin
integer :: kmax
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
      character(len=2 ) :: indfl
!
      if ((kini.eq.1).or.(kini.eq.2)) then
       form='(/,2x,''creation d''''une frontiere rectangle'',/' &
             //'2x,''----------------------------------'',/' &
             //'2x,''numero de la frontiere   : '',11x,i5/'  &
             //'2x,''numero du domaine        : '',11x,i5/'  &
             //'2x,''definition du rectangle  : '',11x,6i5/' &
             //'2x,''type de surface          : '',14x,a/'   &
             //'2x,''cle d''''initialisation     : '',11x,i5)'
      elseif (kini.eq.0) then
       form='(/,2x,''creation d''''une frontiere par'',' &
             //','' fichier d''''indices'',/' &
             //'2x,''----------------------------------------------'',/' &
             //'2x,''fichier d''''indices'',/' &
             //'2x,''numero de la frontiere   : '',11x,i5/'  &
             //'2x,''numero du domaine        : '',11x,i5/'  &
             //'2x,''definition du rectangle  : '',11x,6i5/' &
             //'2x,''type de surface          : '',14x,a/'   &
             //'2x,''cle d''''initialisation     : '',11x,i5)'
      endif
      write(imp,form) mfb,l,imin,imax,jmin,jmax,kmin,kmax,indfl,kini
!
      return
      end
