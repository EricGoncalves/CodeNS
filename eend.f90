module mod_eend
  implicit none
contains
  subroutine eend
    implicit none
!
!***********************************************************************
!
!     ACT
!_A    Action d'arret d'execution du logiciel.
!
!***********************************************************************
!
!$OMP MASTER
    stop
!$OMP END MASTER
!
    return
  end subroutine eend
end module mod_eend
