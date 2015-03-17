module mod_reel
  implicit none
contains
  subroutine reel(mot,imot,r,kerr)
    implicit none
    integer          :: imot,kerr
    double precision :: r
!
!***********************************************************************
!
!     ACT
!_A    Lecture d'un reel dans une variable character.
!
!***********************************************************************
!
    character(len=32) ::  mot
    character(len=7 ) :: formatm
    character(len=2 ) :: longm
!
    write(longm,'(i2)') imot
    formatm='(e'//longm//'.0)'
    kerr=0
    read(mot,formatm,err=100) r
    kerr=1
100 continue
!
    return
  end subroutine reel
end module mod_reel
