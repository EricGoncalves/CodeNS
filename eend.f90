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
    stop
!
    return
  end subroutine eend
end module mod_eend
