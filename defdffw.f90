module mod_defdffw
  implicit none
contains
  subroutine defdffw
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dffw.
!
!***********************************************************************
!
    use kcle
    use definition
    use chainecarac
    implicit none
!
!-----------------------------------------------------------------------
!
!$OMP MASTER
    equat='eu3d '
    kequat=1
!
    klomg=0
    kklomg=1
!
    omg=0.
    komg=1
!$OMP END MASTER
!
    return
  end subroutine defdffw
end module mod_defdffw
