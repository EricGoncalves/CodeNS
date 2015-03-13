module mod_tcmd_dfpmcfg
implicit none
contains
      subroutine tcmd_dfpmcfg(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dfpmcfg.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use kcle
      use chainecarac
      use maillage
      use boundary 
use mod_valenti
implicit none
integer :: imot
integer :: nmot
integer :: icmt
integer :: nm
integer :: no
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
      if(knba.eq.2) knba=3
!
      if(nmot.eq.2)then
        comment=cb
        call synterr(mot,imot,2,comment)
      endif
!
      if(nmot.gt.2) then
        nm=2
        do no=1,mtbx
        nm=nm+1
        if(nmot.lt.nm) then
          comment=ci
          call synterr(mot,imot,nmot,comment)
        else
          call valenti(mot,imot,nm,nba(no),knba)
        endif
        enddo
      else
        comment=cs
        call synterr(mot,imot,nm,comment)
      endif
!
      return
      end subroutine
end module
