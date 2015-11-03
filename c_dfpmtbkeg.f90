module mod_c_dfpmtbkeg
  implicit none
contains
  subroutine c_dfpmtbkeg(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfpmtbkeg.
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
    use mod_b1_dfpmtbkeg
    use mod_tcmd_dfpmtbkeg
    implicit none
    integer          :: imot,nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    dimension imot(nmx)
!
    call tcmd_dfpmtbkeg(mot,imot,nmot)
!
    if(kimp.ge.1) then
       call b1_dfpmtbkeg
    endif
!
    return
  end subroutine c_dfpmtbkeg
end module mod_c_dfpmtbkeg
