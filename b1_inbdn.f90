module mod_b1_inbdn
  implicit none
contains
  subroutine b1_inbdn( &
       lmfb,lmfbd,kibdn)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'inbdn'.
!
!     INP
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    kibdn      : arg int              ; cle initialisat. normales a la front
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    implicit none
    integer          ::       kibdn,      lmfbd,lmfb(lmfbd),       long,     longm1
    integer          ::          nm,      nmult,         nr,      nrest
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=2 ) :: nlg,nlm,nrr,nrm
!
!
    long=6
    longm1=long-1
    nrest=mod(lmfbd,long)
    nmult=(lmfbd-nrest)/long
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
          form='(/,2x,''calcul des normales a une frontiere'',/' &
               //'2x,''-----------------------------------'',/' &
               //'2x,''numero de la frontiere   : '',11x,i5/' &
               //'2x,''cle calcul des normales  : '',11x,i5)'
       else
          form='(/,2x,''donnees de base d''''une frontiere'',/' &
               //'2x,''-------------------------------'',/' &
               //'2x,''numero de la frontiere   : '',11x,i5,' &
               //nrr//'(''       puis'',i5)/' &
               //'2x,''cle calcul des normales  : '',11x,i5)'
       endif
    else if (nm.eq.0) then
       form='(/,2x,''donnees de base d''''une frontiere'',/' &
            //'2x,''-------------------------------'',/' &
            //'2x,''numero de la frontiere   : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //'29x,'//nrr//'(''       puis'',i5)/' &
            //'2x,''cle calcul des normales  : '',11x,i5)'
    else
       form='(/,2x,''donnees de base d''''une frontiere'',/' &
            //'2x,''-------------------------------'',/' &
            //'2x,''numero de la frontiere   : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
            //'29x,'//nrr//'(''       puis'',i5)/' &
            //'2x,''cle calcul des normales  : '',11x,i5)'
    endif
!
    write(imp,form) lmfb,kibdn
!
    return
  end subroutine b1_inbdn
end module mod_b1_inbdn
