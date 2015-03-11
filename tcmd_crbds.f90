module mod_tcmd_crbds
implicit none
contains
      subroutine tcmd_crbds( &
                 mot,imot,nmot, &
                 mfbe,kini,l, &
                 imin,imax,jmin,jmax,kmin,kmax, &
                 indmf)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action crbds.
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
integer :: mfbe
integer :: kini
integer :: l
integer :: imin
integer :: imax
integer :: jmin
integer :: jmax
integer :: kmin
integer :: kmax
integer :: icmt
integer :: im
integer :: kval
integer :: nm
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
      character(len=2 ) :: indmf
      dimension imot(nmx)
!
      do icmt=1,32
      comment(icmt:icmt)=' '
      enddo
      kval=0
!
      nm=3
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,mfbe,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,kini,kval)
      endif
!
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
        call valenti(mot,imot,nm,imin,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,imax,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,jmin,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,jmax,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,kmin,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ci
        call synterr(mot,imot,nmot,comment)
      else
        call valenti(mot,imot,nm,kmax,kval)
      endif
!
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ch
        call synterr(mot,imot,nmot,comment)
      else
        if(imot(nm).gt.2)then
          comment=cc
          call synterr(mot,imot,nm,comment)
        else
          indmf(1:2)='  '
          do im=1,imot(nm)
          indmf(im:im)=mot(nm)(im:im)
          enddo
        endif
      endif
!
      return
      end subroutine
end module
