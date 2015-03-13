module mod_tcmd_cpbd
implicit none
contains
      subroutine tcmd_cpbd(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action cpbd.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use chainecarac
      use kcle
      use boundary
use mod_valenti
implicit none
integer :: imot
integer :: nmot
integer :: icmt
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
      if(kkexl.eq.2) kkexl=3
!
      if(nmot.eq.2)then
        comment=cb
        call synterr(mot,imot,2,comment)
      endif
!
      if(nmot.gt.2) then
        nm=2
        nm=nm+1
        if(nmot.lt.nm) then
          comment=ci
          call synterr(mot,imot,nmot,comment)
        else
          call valenti(mot,imot,nm,kexl,kkexl)
        end if
      else
        comment=cs
        call synterr(mot,imot,nm,comment)
      endif
!
      return
      end subroutine
end module
