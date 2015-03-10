module mod_tcmd_crdms
implicit none
contains
      subroutine tcmd_crdms( &
                 mot,imot,nmot, &
                 l,ni,nj,nk)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action crdms.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
use mod_synterr
use mod_valenti
implicit none
integer :: imot
integer :: nmot
integer :: l
integer :: ni
integer :: nj
integer :: nk
integer :: icmt
integer :: kval
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
      kval=0
!
      nm=3
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,l,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,ni,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,nj,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,nk,kval)
      endif
!
      return
      end
end module
