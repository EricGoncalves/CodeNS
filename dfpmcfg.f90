module mod_dfpmcfg
  implicit none
contains
  subroutine dfpmcfg
!
!***********************************************************************
!
!     ACT
!_A    Definition des parametres de la configuration par la donnee de
!_A    l'ordre d'application des conditions aux limites.
!
!     INP
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use maillage
    use boundary
    implicit none
    integer          :: mfbe,  no
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!$OMP SIMD
    do no=1,mtbx
       mfbe=nba(no)
       nba(no)=nfei(mfbe)
    enddo
!
!$OMP END MASTER
    return
  end subroutine dfpmcfg
end module mod_dfpmcfg
