module mod_b1_inbdb
  implicit none
contains
  subroutine b1_inbdb( &
       lmfb,lmfbd,cl,kibdb, &
       ibdcst,ibdcfl,ibddim,nvbc,vbc)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'inbdb'.
!
!     INP
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    cl         : arg char             ; type de cond lim a appliquer
!_I    kibdb      : arg int              ; cle initialisation de base front
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    implicit none
    integer          ::      ibdcfl,     ibdcst,     ibddim,      kibdb,lmfbd
    integer          :: lmfb(lmfbd),       long,     longm1,         nm,      nmult
    integer          ::         npd,         nr,      nrest,       nvbc
    double precision :: vbc(ista*lsta)
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=8 ) :: cadim
    character(len=4 ) :: cl
    character(len=2 ) :: nlg,nlm,nrr,nrm
!
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
          form='(/,2x,''donnees de base d''''une frontiere'',/' &
               //'2x,''-------------------------------'',/' &
               //'2x,''numero de la frontiere   : '',11x,i5/' &
               //'2x,''cond limite a appliquer  : '',12x,a/' &
               //'2x,''cle calcul de no cell int: '',11x,i5)'
       else
          form='(/,2x,''donnees de base d''''une frontiere'',/' &
               //'2x,''-------------------------------'',/' &
               //'2x,''numero de la frontiere   : '',11x,i5,' &
               //nrr//'(''       puis'',i5)/' &
               //'2x,''cond limite a appliquer  : '',12x,a/' &
               //'2x,''cle calcul de no cell int: '',11x,i5)'
       endif
    else if (nm.eq.0) then
       form='(/,2x,''donnees de base d''''une frontiere'',/' &
            //'2x,''-------------------------------'',/' &
            //'2x,''numero de la frontiere   : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //'29x,'//nrr//'(''       puis'',i5)/' &
            //'2x,''cond limite a appliquer  : '',12x,a/' &
            //'2x,''cle calcul de no cell int: '',11x,i5)'
    else
       form='(/,2x,''donnees de base d''''une frontiere'',/' &
            //'2x,''-------------------------------'',/' &
            //'2x,''numero de la frontiere   : '',11x,i5,' &
            //nlm//'(''       puis'',i5),/' &
            //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
            //'29x,'//nrr//'(''       puis'',i5)/' &
            //'2x,''cond limite a appliquer  : '',12x,a/' &
            //'2x,''cle calcul de no cell int: '',11x,i5)'
    endif
!
    write(imp,form) lmfb,cl,kibdb
!
    if (ibdcst.ne.0) then
       form='(2x,''numero d''''etat            : '',11x,i5)'
       write(imp,form) ibdcst
    endif
!
    if (ibdcfl.ne.0) then
       form='(2x,''no de reference sur fic  : '',11x,i5)'
       write(imp,form) ibdcfl
    endif
!
    if (nvbc.ne.0) then
       if (ibddim.eq.0) cadim='sur ref '
       if (ibddim.eq.1) cadim='sur usi '
!
       long=3
       longm1=long-1
       nrest=mod(nvbc,long)
       nmult=(nvbc-nrest)/long
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
             form='(2x,''donnee physique '',a,'' : '',4x,e12.5)'
          else
             form='(2x,''donnee physique '',a,'' : '',4x,e12.5,' &
                  //nrr//'(''       puis'',e12.5))'
          endif
       else if (nm.eq.0) then
          form='(2x,''donnee physique '',a,'' : '',4x,e12.5,' &
               //nlm//'(''       puis'',e12.5),/' &
               //'29x,'//nrr//'(''       puis'',e12.5))'
       else
          form='(2x,''donnee physique '',a,'' : '',4x,e12.5,' &
               //nlm//'(''       puis'',e12.5),/' &
               //nrm//'(29x,'//nlg//'(''       puis'',e12.5)/)' &
               //'29x,'//nrr//'(''       puis'',e12.5))'
       endif
!
       write(imp,form) cadim,(vbc(npd),npd=1,nvbc)
    endif
!
    return
  end subroutine b1_inbdb
end module mod_b1_inbdb
