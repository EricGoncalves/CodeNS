module mod_b1_secpfw
  implicit none
contains
  subroutine b1_secpfw(lgr,lgrd)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'secpfw'.
!
!     INP
!_I    imp        : com int    ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use kcle
    use schemanum
    use maillage
    use mod_convich
    implicit none
  integer          ::    kncycl,     lgrd,lgr(lgrd),     long,   longm1
  integer          ::        nm,    nmult,       nr,    nrest
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=24) ::  ckvn,cncyresi,cncysave,cncyexpl,cdiscsv,cncycl
    character(len=2 ) :: nlg,nlm,nrr,nrm
!
    call convich(kkvn,ckvn)
    call convich(kncyresi,cncyresi)
    call convich(kncysave,cncysave)
    call convich(kncyexpl,cncyexpl)
    call convich(kdiscsv,cdiscsv)
    kncycl=kncycle(lgr(1))
    call convich(kncycl,cncycl)
!
    form='(/,2x,''environnement de la commande "compute flow"'',/' &
         //'2x,''-------------------------------------------'',/' &
         //'2x,''verification de metrique : '',11x,i5,2x,a/' &
         //'2x,''frequence residu         : '', 8x,i8,2x,a/' &
         //'2x,''frequence exploitation ut: '', 8x,i8,2x,a/' &
         //'2x,''frequence sauvegarde     : '', 8x,i8,2x,a/' &
         //'2x,''discretisation sauvegarde: '',12x,a,2x,a)'
    write(imp,form) kvn,ckvn, &
         ncyresi,cncyresi, &
         ncyexpl,cncyexpl, &
         ncysave,cncysave, &
         discsv,cdiscsv
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
          form='(2x,''realisation'',/' &
               //'2x,''niveau de grille         : '',11x,i5/' &
               //'2x,''ncycle                   : '',11x,i5,2x,a)'
       else
          form='(2x,''realisation'',/' &
               //'2x,''niveau de grille         : '',11x,i5,' &
               //nrr//'(''       puis'',i5)/' &
               //'2x,''ncycle                   : '',11x,i0,2x,a)'
       endif
    else if (nm.eq.0) then
       form='(2x,''realisation'',/' &
            //'2x,''niveau de grille         : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //'29x,'//nrr//'(''       puis'',i5)/' &
            //'2x,''ncycle                   : '',11x,i0,2x,a)'
    else
       form='(2x,''realisation'',/' &
            //'2x,''niveau de grille         : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
            //'29x,'//nrr//'(''       puis'',i5)/' &
            //'2x,''ncycle                   : '',11x,i0,2x,a)'
    endif
!
    write(imp,form) lgr,ncycle(lgr(1)),cncycl
!
    return
  end subroutine b1_secpfw
end module mod_b1_secpfw
