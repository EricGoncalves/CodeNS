module mod_defcpbd
  implicit none
contains
  subroutine defcpbd
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    cpbd.
!
!     COM
!_C    Par defaut on applique les conditions aux limites physiques avant
!_C    l'exploitation des resultats.
!
!***********************************************************************
!
    use kcle
    use boundary
    implicit none
!
!-----------------------------------------------------------------------
!
!$OMP MASTER
    kexl=1
    kkexl=1
!$OMP END MASTER
!
    return
  end subroutine defcpbd
end module mod_defcpbd
