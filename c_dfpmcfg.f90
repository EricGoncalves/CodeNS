module mod_c_dfpmcfg
  implicit none
contains
  subroutine c_dfpmcfg(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfpmcfg.
!
!     INP
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                                   en fct du numero externe
!
!     I/O
!_/    nba        : com int (mtb       ) ; rang de traitement d'une front
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use boundary
    use maillage
    use mod_b1_dfpmcfg

    use mod_dfpmcfg

    use mod_tcmd_dfpmcfg
    implicit none
    integer          :: imot,nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    dimension imot(nmx)
!
    call tcmd_dfpmcfg(mot,imot,nmot)
!
    if(kimp.ge.1) then
       call b1_dfpmcfg
    endif
!
    call dfpmcfg
!
    return
  end subroutine c_dfpmcfg
end module mod_c_dfpmcfg
