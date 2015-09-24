module mod_b1_dfpmimd
  implicit none
contains
  subroutine b1_dfpmimd( &
       ldom,ldomd, &
       lgr,lgrd)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfpmimd'.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    kmf        : arg int              ; cle phase implicite
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use schemanum
    use kcle
    use chainecarac
    use mod_convich
    implicit none
    integer          ::       ldomd,       lgrd,ldom(ldomd),  lgr(lgrd),         lm
    integer          ::        long,     longm1,         nm,      nmult,         nr
    integer          ::       nrest
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=24) ::  ckmf,clmax
    character(len=2 ) :: nlg,nlm,nrr,nrm
!
!
    lm=ldom(1)+(lgr(1)-1)*lz
!
    call convich(kkmf(lm),ckmf)
    call convich(klmax(lm),clmax)
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
          form='(/,2x,''parametres de resolution implicite'',/' &
               //'2x,''----------------------------------'',/' &
               //'2x,''numero du domaine        : '',11x,i5)'
       else
          form='(/,2x,''parametres de resolution implicite'',/' &
               //'2x,''----------------------------------'',/' &
               //'2x,''numero du domaine        : '',11x,i5,' &
               //nrr//'(''       puis'',i5))'
       endif
    else if (nm.eq.0) then
       form='(/,2x,''parametres de resolution implicite'',/' &
            //'2x,''----------------------------------'',/' &
            //'2x,''numero du domaine        : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //'29x,'//nrr//'(''       puis'',i5))'
    else
       form='(/,2x,''parametres de resolution implicite'',/' &
            //'2x,''----------------------------------'',/' &
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
          form='(2x,''niveau de grille         : '',11x,i5/' &
               //'2x,''kmf                     : '',11x,i5,2x,a)'
       else
          form='(2x,''niveau de grille         : '',11x,i5,' &
               //nrr//'(''       puis'',i5)/' &
               //'2x,''kmf                     : '',11x,i5,2x,a)'
       endif
    else if (nm.eq.0) then
       form='(2x,''niveau de grille         : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //'29x,'//nrr//'(''       puis'',i5)/' &
            //'2x,''kmf                     : '',11x,i5,2x,a)'
    else
       form='(2x,''niveau de grille         : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
            //'29x,'//nrr//'(''       puis'',i5)/' &
            //'2x,''kmf                     : '',11x,i5,2x,a)'
    endif
!
    write(imp,form) lgr,kmf(lm),ckmf
!
    return
  end subroutine b1_dfpmimd
end module mod_b1_dfpmimd
