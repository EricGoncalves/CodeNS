module mod_b1_dfpmcfg
  implicit none
contains
  subroutine b1_dfpmcfg
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfpmcfg'.
!
!     INP
!_I    mtbx       : arg int              ; nbr total de frontieres
!_I    nba        : arg int (mtb       ) ; rang de traitement d'une front
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use maillage
    use boundary
    use kcle
    use mod_convich
    implicit none
    integer          ::   long,longm1,    nm, nmult,    no
    integer          ::     nr, nrest
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=24) ::  cnba
    character(len=2 ) :: nlg,nlm,nrr,nrm
!
    call convich(knba,cnba)
!
    long=6
    longm1=long-1
    nrest=mod(mtbx,long)
    nmult=(mtbx-nrest)/long
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
          form='(/,2x,''parametres de configuration'',/' &
               //'2x,''---------------------------'',/' &
               //'2x,''fr ds ordre traitement   : '',11x,i5)'
       else
          form='(/,2x,''parametres de configuration'',/' &
               //'2x,''---------------------------'',/' &
               //'2x,''fr ds ordre traitement   : '',11x,i5,' &
               //nrr//'(''       puis'',i5))'
       endif
    else if (nm.eq.0) then
       form='(/,2x,''parametres de configuration'',/' &
            //'2x,''---------------------------'',/' &
            //'2x,''fr ds ordre traitement   : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //'29x,'//nrr//'(''       puis'',i5))'
    else
       form='(/,2x,''parametres de configuration'',/' &
            //'2x,''---------------------------'',/' &
            //'2x,''fr ds ordre traitement   : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
            //'29x,'//nrr//'(''       puis'',i5))'
    endif
!
    write(imp,form) (nba(no),no=1,mtbx)
!
    return
  end subroutine b1_dfpmcfg
end module mod_b1_dfpmcfg
