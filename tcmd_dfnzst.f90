      subroutine tcmd_dfnzst( &
                 mot,imot,nmot, &
                 nonzst)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dfnzst.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
!
!
      do icmt=1,32
       comment(icmt:icmt)=' '
      enddo
!
      kval=0
!
      if(nmot.eq.2)then
        comment=cb
        call synterr(mot,imot,2,comment)
      endif
!
      if(nmot.gt.2) then
        nm=2
!
        nm=nm+1
          call valenti(mot,imot,nm,nonzst,kval)
      endif
!
      return
      end
