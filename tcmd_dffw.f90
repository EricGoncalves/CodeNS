      subroutine tcmd_dffw(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dffw.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use kcle
      use definition
implicit none
integer :: imot
integer :: nmot
integer :: icmt
integer :: im
integer :: nm
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
!
      do icmt=1,32
      comment(icmt:icmt)=' '
      enddo
!
      if(kequat.eq.2) kequat=3
      if(kklomg.eq.2) kklomg=3
      if(komg.eq.2) komg=3
!
      if(nmot.eq.2)then
        comment=cb
        call synterr(mot,imot,2,comment)
      endif
!
      if(nmot.gt.2) then
        nm=2
      do while(nm.lt.nmot)
        nm=nm+1
        if((imot(nm).eq.5).and.(mot(nm).eq.'equat')) then
          nm=nm+1
          if(imot(nm).gt.7)then
            comment=cc
            call synterr(mot,imot,nm,comment)
          else
            equat(1:7)='       '
            do im=1,imot(nm)
            equat(im:im)=mot(nm)(im:im)
            enddo
            kequat=2
          endif
        else if((imot(nm).eq.9).and.(mot(nm).eq.'krotation')) then
          nm=nm+1
          call valenti(mot,imot,nm,klomg,kklomg)
        elseif((imot(nm).eq.9).and.(mot(nm).eq.'vrotation')) then
          nm=nm+1
          call valreel(mot,imot,nm,omg,komg)
        elseif(imot(nm).eq.0) then
          comment=cs
          call synterr(mot,imot,nm,comment)
        else
          comment=cb
          call synterr(mot,imot,nm,comment)
        endif
       enddo
      endif
!
      return
      end
