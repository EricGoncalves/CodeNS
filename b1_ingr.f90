module mod_b1_ingr
implicit none
contains
      subroutine b1_ingr( &
                 ldom,ldomd,king)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'ingr'.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    king       : arg int              ; cle initialisation maillage
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!
      use sortiefichier
implicit none
integer :: ldom
integer :: ldomd
integer :: king
integer :: long
integer :: longm1
integer :: nm
integer :: nmult
integer :: nr
integer :: nrest
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
      character(len=2 ) :: nlg,nlm,nrr,nrm
      dimension ldom(ldomd)
!
      long=6
      longm1=long-1
      nrest=mod(ldomd,long)
      nmult=(ldomd-nrest)/long
      if(nrest.eq.0) then
        nr=long
        nm=nmult-2
      else
        nr=nrest
        nm=nmult-1
      endif
      if (nm.eq.-1) nr=nr-1
!
      write(nlg,'(i2)') long
      write(nlm,'(i2)') longm1
      write(nrr,'(i2)') nr
      write(nrm,'(i2)') nm
!
      if (nm.eq.-1) then
        if(nr.eq.0) then
          form='(/,2x,''initialisation de maillage'',/' &
                //'2x,''--------------------------'',/' &
                //'2x,''numero du domaine           : '',11x,i5/' &
                //'2x,''cle d''''initialisation     : '',11x,i5)'
        else
          form='(/,2x,''initialisation de maillage'',/' &
                //'2x,''--------------------------'',/' &
                //'2x,''numero du domaine           : '',11x,i5,' &
                //nrr//'(''       puis'',i5)/' &
                //'2x,''cle d''''initialisation     : '',11x,i5)'
        endif
      elseif (nm.eq.0) then
        form='(/,2x,''initialisation de maillage'',/' &
           //'2x,''--------------------------'',/' &
           //'2x,''numero du domaine           : '',11x,i5,' &
           //nlm//'(''       puis'',i5),/' &
           //'29x,'//nrr//'(''       puis'',i5)/' &
           //'2x,''cle d''''initialisation     : '',11x,i5)'
      else
        form='(/,2x,''initialisation de maillage'',/' &
           //'2x,''--------------------------'',/' &
           //'2x,''numero du domaine           : '',11x,i5,' &
           //nlm//'(''       puis'',i5),/' &
           //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
           //'29x,'//nrr//'(''       puis'',i5)/' &
           //'2x,''cle d''''initialisation     : '',11x,i5)'
      endif
!
      write(imp,form) ldom,king
!
      return
      end subroutine
end module
