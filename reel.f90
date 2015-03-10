      subroutine reel(mot,imot,r,kerr)
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
      end
