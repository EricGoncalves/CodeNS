module mod_c_dftl1
  implicit none
contains
  subroutine c_dftl1(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dftl1.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use mod_tcmd_dftl1
    implicit none
    integer          :: imot(nmx),     nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
!
    call tcmd_dftl1(mot,imot,nmot)
!
!      if (kimp.ge.1) then
!            call b1_dftl1
!      endif
!
    return
  end subroutine c_dftl1
end module mod_c_dftl1
