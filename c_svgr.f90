module mod_c_svgr
  implicit none
contains
  subroutine c_svgr( &
       mot,imot,nmot, &
       l,x,y,z, &
       tn1,tn2,tn3)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action svgr.
!
!_I    l          : arg int              ; numero de domaine
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_L    tn1        : arg real(ip00      ) ; tableau de travail
!_L    tn2        : arg real(ip00      ) ; tableau de travail
!_L    tn3        : arg real(ip00      ) ; tableau de travail
!_I    disc       : arg char             ; changement de discretisation (centre
!_I                                        /noeud)
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!
!     COM
!_C    Les donnees suvegardees ici sont celles qui peuvent etre demandees
!_C    pour l' initialisation de donnees d'un domaine par ingr.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use mod_svgr
    use mod_tcmd_svgr
    implicit none
    integer          :: imot(nmx),        l,     nmot
    double precision :: tn1(ip00),tn2(ip00),tn3(ip00),  x(ip21),  y(ip21)
    double precision ::   z(ip21)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    character(len=4 ) :: disc
!$OMP MASTER
!
!
    call tcmd_svgr( &
         mot,imot,nmot, &
         disc)
!
!      if ((kimp.ge.1).and.(l.eq.1)) then
!            call b1_svgr(disc)
!      endif
!
    call svgr( &
         l,x,y,z, &
         tn1,tn2,tn3, &
         disc)
!
!$OMP END MASTER
    return
  end subroutine c_svgr
end module mod_c_svgr
