module mod_c_dfnm
  implicit none
contains
  subroutine c_dfnm(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfnm.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use mod_tcmd_dfnm
    use mod_b1_dfnm

    implicit none
    integer          :: imot(nmx),     nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
!
    call tcmd_dfnm(mot,imot,nmot)
!
    if (kimp.ge.1) then
       call b1_dfnm
    endif
!
    return
  end subroutine c_dfnm
end module mod_c_dfnm
