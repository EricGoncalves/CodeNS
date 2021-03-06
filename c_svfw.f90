module mod_c_svfw
  implicit none
contains
  subroutine c_svfw( &
       mot,imot,nmot, &
       l,v,mut,utau, &
       ncin,ncbd, &
       tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action svfw.
!
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    disc       : arg char             ; changement de discretisation (centre
!_I                                        /noeud)
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_L    tn1        : arg real(ip00      ) ; tableau de travail
!_L    tn2        : arg real(ip00      ) ; tableau de travail
!_L    tn3        : arg real(ip00      ) ; tableau de travail
!_L    tn4        : arg real(ip00      ) ; tableau de travail
!_L    tn5        : arg real(ip00      ) ; tableau de travail
!_L    tn6        : arg real(ip00      ) ; tableau de travail
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!
!
!     COM
!_C    Les donnees sauvegardees ici sont celles qui peuvent etre demandees
!_C    pour l' initialisation de donnees d'un domaine par infw.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use mod_svfw
    use mod_tcmd_svfw
    use mod_b1_svfw
    use mod_mpi
    implicit none
    integer          ::  imot(nmx),         l,ncbd(ip41),ncin(ip41),      nmot
    double precision ::    mut(ip12),   tn1(ip00),   tn2(ip00),   tn3(ip00),   tn4(ip00)
    double precision ::    tn5(ip00),   tn6(ip00),   tn7(ip00),   tn8(ip00),  utau(ip42)
    double precision :: v(ip11,ip60)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    character(len=4 ) :: disc
!
!
    call tcmd_svfw( &
         mot,imot,nmot, &
         disc)
!
    if ((kimp.ge.1).and.(bl_to_bg(l).eq.1)) then
       call b1_svfw(disc)
    endif
!
    call svfw( &
         l,v,mut,utau, &
         disc,ncin,ncbd, &
         tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8)
!
    return
  end subroutine c_svfw
end module mod_c_svfw
