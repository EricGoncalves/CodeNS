module mod_defdfgm
  implicit none
contains
  subroutine defdfgm
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfgm.
!
!     COM
!_C    Par defaut configuration sans machine tournante et
!_C    sans periodicite.
!
!***********************************************************************
!
    use kcle
    use chainecarac
    use definition
    implicit none
!
!-----------------------------------------------------------------------
!
!$OMP MASTER
    config='tuy '
    kconfig=1
!
    perio=0.
    kperio=1
!$OMP END MASTER
!
    return
  end subroutine defdfgm
end module mod_defdfgm
