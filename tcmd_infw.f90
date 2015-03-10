      subroutine tcmd_infw( &
                 mot,imot,nmot, &
                 ldom,ldomd,kina)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action infw.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use maillage
      use kcle
      use modeleturb
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
      dimension ldom(nobj)
!
      do icmt=1,32
       comment(icmt:icmt)=' '
      enddo
      kval=0
!
      nm=3
      nm=nm+1
      if(nmot.lt.nm) then
        comment=cm
        call synterr(mot,imot,nmot,comment)
      else
            call vallent(mot,imot,nm,ldom,ldomd,lzx,klzx)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,kina,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
       return
      elseif((imot(nm).eq.6).and.(mot(nm).eq.'keinit')) then
          nm=nm+1
          if(nmot.lt.nm) then
            comment=ci
            call synterr(mot,imot,nmot,comment)
          else
          call valenti(mot,imot,nm,keinit,kval)
          endif
      else
         comment=cs
         call synterr(mot,imot,nm,comment)
      endif
!
      return
      end
