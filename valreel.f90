module mod_valreel
implicit none
contains
      subroutine valreel(mot,imot,nm,rree,krree)
!
!***********************************************************************
!
!     ACT
!_A    Affectation de sa valeur au reel rree.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
use mod_synterr
use mod_reel
implicit none
integer :: imot
integer :: nm
double precision :: rree
integer :: krree
integer :: icmt
integer :: kerr
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
      call reel(mot(nm),imot(nm),rree,kerr)
      if(kerr.eq.0)then
        comment=cr
        call synterr(mot,imot,nm,comment)
      else
        krree=2
      endif
!
      return
      end subroutine
end module
