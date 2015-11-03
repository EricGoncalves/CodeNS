module mod_def
  implicit none
contains
  subroutine def
    use mod_defcpbd
    use mod_defdftl1
    use mod_defdfgm
    use mod_defintn
    use mod_defdfnzst
    use mod_defdfph
    use mod_defdffw
    use mod_defdfst
    use mod_defdfpmimd
    use mod_defsecpfw
    use mod_defdfpmdsd
    use mod_defdfpmdtg
    use mod_defdfpmtbn
    use mod_defdfpmcfg
    use mod_defdfnm
    use mod_defdfpmdtd
    implicit none
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut pour tout le code.
!
!***********************************************************************
!
    call defcpbd
    call defdffw
    call defdfgm
    call defdfnm
    call defdfnzst
    call defdfph
!    call defdfpmcfg
!    call defdfpmdtd
    call defdfpmdtg
!    call defdfpmdsd
!    call defdfpmimd
    call defdfpmtbn
    call defdfst
    call defdftl1
    call defintn
    call defsecpfw
!
    return
  end subroutine def
end module mod_def
