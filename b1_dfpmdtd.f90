      subroutine b1_dfpmdtd( &
                 ldom,ldomd, &
                 lgr,lgrd)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfpmdtd'.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    eta        : arg real             ; nombre de CFL
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use sortiefichier
   use schemanum
   use kcle
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
      character(len=24) ::  ceta
      character(len=2 ) :: nlg,nlm,nrr,nrm
!
      dimension ldom(ldomd)
      dimension lgr(lgrd)
!
      lm=ldom(1)+(lgr(1)-1)*lz
!
      call convich(keta(lm),ceta)
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
          form='(/,2x,''nombre CFL'',/' &
                //'2x,''----------'',/' &
                //'2x,''numero du domaine        : '',11x,i5)'
        else
          form='(/,2x,''nombre CFL'',/' &
                //'2x,''----------'',/' &
                //'2x,''numero du domaine        : '',11x,i5,' &
                //nrr//'(''       puis'',i5))'
        endif
      else if (nm.eq.0) then
          form='(/,2x,''nombre CFL'',/' &
                //'2x,''----------'',/' &
              //'2x,''numero du domaine        : '',11x,i5,' &
              //nlm//'(''       puis'',i5),/' &
              //'29x,'//nrr//'(''       puis'',i5))'
      else
          form='(/,2x,''nombre CFL'',/' &
                //'2x,''----------'',/' &
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
               //'2x,''eta                      : '',4x,e12.5,2x,a)'
        else
          form='(2x,''niveau de grille         : '',11x,i5,' &
               //nrr//'(''       puis'',i5)/' &
               //'2x,''eta                      : '',4x,e12.5,2x,a)'
        endif
      else if (nm.eq.0) then
        form='(2x,''niveau de grille         : '',11x,i5,' &
             //nlm//'(''       puis'',i5),/' &
             //'29x,'//nrr//'(''       puis'',i5)/' &
             //'2x,''eta                      : '',4x,e12.5,2x,a)'
      else
        form='(2x,''niveau de grille         : '',11x,i5,' &
             //nlm//'(''       puis'',i5),/' &
             //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
             //'29x,'//nrr//'(''       puis'',i5)/' &
             //'2x,''eta                      : '',4x,e12.5,2x,a)'
      endif
!
      write(imp,form) lgr,eta(lm),ceta
!
      return
      end
