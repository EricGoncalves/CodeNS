      subroutine entier(mot,imot,i,kerr)
!
!***********************************************************************
!
!     ACT
!_A    Lecture d'un entier dans une variable character.
!
!***********************************************************************
!
      character *32 mot
      character *4 formatm
      character *1 longm
!
       write(longm,'(i1)') imot
       formatm='(i'//longm//')'
       kerr=0
       read(mot,formatm,err=100) i
       kerr=1
  100 continue
!
      return
      end
