      subroutine tcmd_svgr( &
                 mot,imot,nmot, &
                 disc)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action svgr.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
!
!-----------------------------------------------------------------------
!
      character *32 comment
      character *32 mot(nmx)
      character *4 disc
      dimension imot(nmx)
!
      do icmt=1,32
      comment(icmt:icmt)=' '
      enddo
!
      nm=3
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ch
        call synterr(mot,imot,nmot,comment)
      else
        if(imot(nm).gt.4)then
          comment=cc
          call synterr(mot,imot,nm,comment)
        else
          disc(1:4)='    '
          do im=1,imot(nm)
          disc(im:im)=mot(nm)(im:im)
          enddo
        endif
      endif
!
      return
      end
