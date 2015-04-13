module mod_c_dfpmdtg
  implicit none
contains
  subroutine c_dfpmdtg(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfpmdtg.
!
!     INP
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    icytur0    : com int              ; nbr de cycl en deb de calcul au cours
!_I                                        desquelles mut n'est pas mis a jour
!_I    ncyturb    : com int              ; freq en it de mise a jour de mut
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use mod_b1_dfpmdtg

    use mod_tcmd_dfpmdtg
    implicit none
    integer          :: imot(nmx),     nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
!$OMP MASTER
!
    call tcmd_dfpmdtg(mot,imot,nmot)
!
    if(kimp.ge.1) then
       call b1_dfpmdtg
    endif
!
!$OMP END MASTER
    return
  end subroutine c_dfpmdtg
end module mod_c_dfpmdtg
