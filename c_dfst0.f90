module mod_c_dfst0
  implicit none
contains
  subroutine c_dfst0(roam,aam,tam)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfst0.
!
!     INP
!_I    roam       : arg real             ; etat de reference utilisateur
!_I                                        dimensionne, masse volumique d'arret
!_I    aam        : arg real             ; etat de reference utilisateur
!_I                                        dimensionne, vitesse du son d'arret
!_I    tam        : arg real             ; etat de reference utilisateur
!_I                                        dimensionne, temperature d'arret
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!
!     COM
!_C    Cet etat 'am', une fois normalise et note 'a1', est utilise pour
!_C    les conditions d'injection (clidd,clidi,cliti,clivi,clivd) et
!_C    la condition d'adherence (clpari).
!
!***********************************************************************
!
    use sortiefichier
    use mod_b1_dfst0
    implicit none
    double precision ::  aam,roam, tam
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
    if (kimp.ge.1) then
       call b1_dfst0(roam,aam,tam)
    endif
!
!$OMP END MASTER
    return
  end subroutine c_dfst0
end module mod_c_dfst0
