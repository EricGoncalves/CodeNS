module mod_lentier
implicit none
contains
      subroutine lentier(mot,imot,li,lid,lix,klix,kerr)
!
!***********************************************************************
!
!     ACT
!_A    Lecture d'un entier dans une variable chaine caractere.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
use mod_entier
implicit none
integer :: imot
integer :: li
integer :: lid
integer :: lix
integer :: klix
integer :: kerr
integer :: l
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot
      dimension li(nobj)
!
      if ((imot.eq.1).and.(mot.eq.'*')) then
        if(klix.eq.0) then
          kerr=-1
        else
          do l=1,lix
          li(l)=l
          enddo
          lid=lix
          kerr=1
        endif
      else
        call entier(mot,imot,li(1),kerr)
        lid=1
      endif
!
      return
!
      end
end module
