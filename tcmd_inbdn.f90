      subroutine tcmd_inbdn( &
                 mot,imot,nmot, &
                 lmfb,lmfbd,kibdn)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action inbdn.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use maillage
      use kcle
implicit none
integer :: imot
integer :: nmot
integer :: lmfb
integer :: lmfbd
integer :: kibdn
integer :: icmt
integer :: kval
integer :: nm
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
      dimension lmfb(nobj)
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
        call vallent(mot,imot,nm,lmfb,lmfbd,mtbx,kmtbx)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,kibdn,kval)
      endif
!
      return
      end
