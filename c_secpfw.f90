module mod_c_secpfw
  implicit none
contains
  subroutine c_secpfw(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action secpfw.
!
!***********************************************************************
!
    use para_fige
    use sortiefichier
    use mod_tcmd_secpfw
    use mod_b1_secpfw
    implicit none
  integer          :: imot(nmx),lgr(nobj),     lgrd,     nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
!
!
    call tcmd_secpfw( &
         mot,imot,nmot, &
         lgr,lgrd)
!
    if(kimp.ge.1) then
       call b1_secpfw(lgr,lgrd)
    endif
!
    return
  end subroutine c_secpfw
end module mod_c_secpfw
