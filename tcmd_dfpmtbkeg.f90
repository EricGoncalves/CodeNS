module mod_tcmd_dfpmtbkeg
implicit none
contains
      subroutine tcmd_dfpmtbkeg(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dfpmtbkeg.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use modeleturb
      use schemanum
use mod_valreel
implicit none
integer :: imot
integer :: nmot
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
!
      if(nmot.eq.4)then
        comment=cb
        call synterr(mot,imot,4,comment)
      endif
!
        nm=4
!---
        nm=nm+1
        if((imot(nm).eq.6).and.(mot(nm).eq.'rokinf')) then
          nm=nm+1
                if(nmot.lt.nm) then
                  comment=cr
                  call synterr(mot,imot,nmot,comment)
                else
            call valreel(mot,imot,nm,rokinf,kval)
                endif
      else
         comment=cs
         call synterr(mot,imot,nm,comment)
      endif
!---
        nm=nm+1
        if((imot(nm).eq.6).and.(mot(nm).eq.'roeinf')) then
          nm=nm+1
                if(nmot.lt.nm) then
                  comment=cr
                  call synterr(mot,imot,nmot,comment)
                else
            call valreel(mot,imot,nm,roeinf,kval)
                endif
      else
         comment=cs
         call synterr(mot,imot,nm,comment)
      endif
!---
        nm=nm+1
        if((imot(nm).eq.4).and.(mot(nm).eq.'epsk')) then
          nm=nm+1
                if(nmot.lt.nm) then
                  comment=cr
                  call synterr(mot,imot,nmot,comment)
                else
            call valreel(mot,imot,nm,epsk,kval)
                endif
      else
         comment=cs
         call synterr(mot,imot,nm,comment)
      endif
!---
        nm=nm+1
        if((imot(nm).eq.4).and.(mot(nm).eq.'epse')) then
          nm=nm+1
                if(nmot.lt.nm) then
                  comment=cr
                  call synterr(mot,imot,nmot,comment)
                else
            call valreel(mot,imot,nm,epse,kval)
                endif
      else
         comment=cs
         call synterr(mot,imot,nm,comment)
      endif
!---
        nm=nm+1
        if((imot(nm).eq.5).and.(mot(nm).eq.'rki2t')) then
          nm=nm+1
                if(nmot.lt.nm) then
                  comment=cr
                  call synterr(mot,imot,nmot,comment)
                else
            call valreel(mot,imot,nm,rki2t,kval)
                endif
      else
         comment=cs
         call synterr(mot,imot,nm,comment)
      endif
!---
        nm=nm+1
        if((imot(nm).eq.5).and.(mot(nm).eq.'rki4t')) then
          nm=nm+1
                if(nmot.lt.nm) then
                  comment=cr
                  call synterr(mot,imot,nmot,comment)
                else
            call valreel(mot,imot,nm,rki4t,kval)
                endif
      else
         comment=cs
         call synterr(mot,imot,nm,comment)
      endif
!
      return
      end subroutine
end module
