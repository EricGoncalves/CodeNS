module mod_b1_infw
implicit none
contains
      subroutine b1_infw( &
                 ldom,ldomd,kina)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'infw'.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    kina       : arg int              ; cle initialisation var de calcul
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!
      use sortiefichier
   use chainecarac
   use modeleturb
implicit none
integer :: ldom
integer :: ldomd
integer :: kina
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
      if (nrest.eq.0) then
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
          form='(/,2x,''initialisation aerodynamique'',/' &
                //'2x,''----------------------------'',/' &
                //'2x,''numero du domaine        : '',11x,i5/' &
                //'2x,''cle d''''initialisation     : '',11x,i5)'
        else
          form='(/,2x,''initialisation aerodynamique'',/' &
                //'2x,''----------------------------'',/' &
                //'2x,''numero du domaine           : '',11x,i5,' &
                //nrr//'(''       puis'',i5)/' &
                //'2x,''cle d''''initialisation     : '',11x,i5)'
        endif
      else if (nm.eq.0) then
        form='(/,2x,''initialisation aerodynamique'',/' &
              //'2x,''----------------------------'',/' &
              //'2x,''numero du domaine           : '',11x,i5,' &
              //nlm//'(''       puis'',i5),/' &
              //'29x,'//nrr//'(''       puis'',i5)/' &
              //'2x,''cle d''''initialisation     : '',11x,i5)'
      else
        form='(/,2x,''initialisation aerodynamique'',/' &
              //'2x,''----------------------------'',/' &
              //'2x,''numero du domaine           : '',11x,i5,' &
              //nlm//'(''       puis'',i5),/' &
              //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
              //'29x,'//nrr//'(''       puis'',i5)/' &
              //'2x,''cle d''''initialisation     : '',11x,i5)'
      endif
!
      write(imp,form) ldom,kina
!
      if(equat(6:7).eq.'ke') then
          form='(2x,''type initialisation k-eps: '',11x,i5)'
      write(imp,form) keinit
      endif
!
      return
      end
end module
