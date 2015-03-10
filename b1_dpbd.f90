module mod_b1_dpbd
implicit none
contains
      subroutine b1_dpbd( &
                 lmfb,lmfbd, &
                 lgr,lgrd, &
                 typdat)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dpbd'.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use sortiefichier
implicit none
integer :: lmfb
integer :: lmfbd
integer :: lgr
integer :: lgrd
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
      character(len=32) ::  typdat
      character(len=2 ) :: nlg,nlm,nrr,nrm
!
      dimension lmfb(lmfbd)
      dimension lgr(lgrd)
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
          form='(/,2x,''impression frontiere : '',a/' &
                //'2x,''--------------------'',/' &
                //'2x,''numero de frontiere      : '',11x,i5)'
        else
          form='(/,2x,''impression frontiere : '',a/' &
                //'2x,''--------------------'',/' &
                //'2x,''numero de frontiere      : '',11x,i5,' &
                //nrr//'(''       puis'',i5))'
        endif
      else if (nm.eq.0) then
          form='(/,2x,''impression frontiere : '',a/' &
                //'2x,''--------------------'',/' &
              //'2x,''numero de frontiere      : '',11x,i5,' &
              //nlm//'(''       puis'',i5),/' &
              //'29x,'//nrr//'(''       puis'',i5))'
      else
          form='(/,2x,''impression frontiere : '',a/' &
                //'2x,''--------------------'',/' &
              //'2x,''numero de frontiere      : '',11x,i5,' &
              //nlm//'(''       puis'',i5),/' &
              //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
              //'29x,'//nrr//'(''       puis'',i5))'
      endif
!
      write(imp,form) typdat,lmfb
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
          form='(2x,''niveau de grille         : '',11x,i5)'
        else
          form='(2x,''niveau de grille         : '',11x,i5,' &
               //nrr//'(''       puis'',i5))'
        endif
      else if (nm.eq.0) then
        form='(2x,''niveau de grille         : '',11x,i5,' &
             //nlm//'(''       puis'',i5),/' &
             //'29x,'//nrr//'(''       puis'',i5))'
      else
        form='(2x,''niveau de grille         : '',11x,i5,' &
             //nlm//'(''       puis'',i5),/' &
             //nrm//'(29x,'//nlg//'(''       puis'',i5)/)' &
             //'29x,'//nrr//'(''       puis'',i5))'
      endif
!
      write(imp,form) lgr
!
      return
      end
end module
