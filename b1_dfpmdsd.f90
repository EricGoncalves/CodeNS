module mod_b1_dfpmdsd
implicit none
contains
      subroutine b1_dfpmdsd( &
                 ldom,ldomd, &
                 lgr,lgrd)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfpmdsd'.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    ki2        : arg real             ; coef de dissipation ordre 2
!_I    ki4        : arg real             ; coef de dissipation ordre 4
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use sortiefichier
   use schemanum
   use kcle
use mod_convich
implicit none
integer :: ldom
integer :: ldomd
integer :: lgr
integer :: lgrd
integer :: lm
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
      character(len=24) ::  cki2,cki4
      character(len=2 ) :: nlg,nlm,nrr,nrm
!
      dimension ldom(ldomd)
      dimension lgr(lgrd)
!
      lm=ldom(1)+(lgr(1)-1)*lz
!
      call convich(kki2(lm),cki2)
      call convich(kki4(lm),cki4)
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
          form='(/,2x,''parametres de dissipation'',/' &
                //'2x,''-------------------------'',/' &
                //'2x,''numero du domaine        : '',11x,i5)'
        else
          form='(/,2x,''parametres de dissipation'',/' &
                //'2x,''-------------------------'',/' &
                //'2x,''numero du domaine        : '',11x,i5,' &
                //nrr//'(''       puis'',i5))'
        endif
      else if (nm.eq.0) then
          form='(/,2x,''parametres de dissipation'',/' &
                //'2x,''-------------------------'',/' &
              //'2x,''numero du domaine        : '',11x,i5,' &
              //nlm//'(''       puis'',i5),/' &
              //'29x,'//nrr//'(''       puis'',i5))'
      else
          form='(/,2x,''parametres de dissipation'',/' &
                //'2x,''-------------------------'',/' &
              //'2x,''numero du domaine        : '',11x,i5,' &
              //nlm//'(''       puis'',i5),/' &
              //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
              //'29x,'//nrr//'(''       puis'',i5))'
      endif
!
      write(imp,form) ldom
!
      long=6
      longm1=long-1
      nrest=mod(lgrd,long)
      nmult=(lgrd-nrest)/long
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
          form='(2x,''niveau de grille          : '',11x,i5/' &
               //'2x,''ki2                      : '',4x,e12.5,2x,a/' &
               //'2x,''ki4                      : '',4x,e12.5,2x,a)'
        else
          form='(2x,''niveau de grille         : '',11x,i5,' &
               //nrr//'(''       puis'',i5)/' &
               //'2x,''ki2                      : '',4x,e12.5,2x,a/' &
               //'2x,''ki4                      : '',4x,e12.5,2x,a)'
        endif
      else if (nm.eq.0) then
        form='(2x,''niveau de grille         : '',11x,i5,' &
             //nlm//'(''       puis'',i5),/' &
             //'29x,'//nrr//'(''       puis'',i5)/' &
             //'2x,''ki2                      : '',4x,e12.5,2x,a/' &
             //'2x,''ki4                      : '',4x,e12.5,2x,a)'
      else
        form='(2x,''niveau de grille         : '',11x,i5,' &
             //nlm//'(''       puis'',i5),/' &
             //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
             //'29x,'//nrr//'(''       puis'',i5)/' &
             //'2x,''ki2                      : '',4x,e12.5,2x,a/' &
             //'2x,''ki4                      : '',4x,e12.5,2x,a)'
      endif
!
      write(imp,form) lgr,ki2(lm),cki2, &
                          ki4(lm),cki4
!
      return
      end
end module
