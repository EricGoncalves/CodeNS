module mod_c_nzst
  implicit none
contains
  subroutine c_nzst(roam,aam,tam)
!
!***********************************************************************
!
!     ACT
!_A    Normalisation d'un etat note 'am' par l'etat note 'nz', l'etat
!_A    adimensionne est note 'a1': roa1,aa1,ta1,pa1,ha1.
!
!     INP
!_I    c_nzst     roam       : arg real             ; etat de reference utilisateur
!_I                                                   dimensionne, masse volumique d'arret
!_I    c_nzst     aam        : arg real             ; etat de reference utilisateur
!_I                                                   dimensionne, vitesse du son d'arret
!_I    c_nzst     tam        : arg real             ; etat de reference utilisateur
!_I                                                   dimensionne, temperature d'arret
!_I    c_nzst     kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    c_nzst     gam        : com real             ; rapport des chaleurs specifiques
!_I    c_nzst     rgp        : com real             ; constante des gaz parfaits adim
!_I    c_nzst     tnz        : com real             ; etat pour adimensionnement,
!_I                                                   temperature
!_I    c_nzst     ronz       : com real             ; etat pour adimensionnement,
!_I                                                   masse volumique
!_I    c_nzst     anz        : com real             ; etat pour adimensionnement,
!_I                                                   vitesse du son d'arret
!
!     OUT
!_O    c_nzst     roa1       : arg real             ; etat de reference utilisateur
!_O                                                   adimensionne, masse volumique d'arret
!_O    c_nzst     aa1        : arg real             ; etat de reference utilisateur
!_O                                                   adimensionne, vitesse du son d'arret
!_O    c_nzst     ta1        : arg real             ; etat de reference utilisateur
!_O                                                   adimensionne, temperature d'arret
!_O    c_nzst     pa1        : arg real             ; pression d'arret de l'etat
!_O                                                   de reference utilisateur adimensionne
!_O    c_nzst     ha1        : arg real             ; enthalpie d'arret de l'etat
!_O                                                   de reference utilisateur adimensionne
!
!     COM
!_C    Cet etat 'a1' est utilise pour
!_C    les conditions d'injection (clidd,clidi,cliti,clivi,clivd) et
!_C    la condition d'adherence (clpari).
!
!***********************************************************************
!
    use sortiefichier
    use proprieteflu
    use definition
    implicit none
    double precision ::  aam,roam, tam
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
    roa1=roam/ronz
    aa1=aam/anz
    ta1=tam/tnz
    pa1=roa1*ta1/gam
    ha1=cp*ta1
!
!$OMP END MASTER
    return
  end subroutine c_nzst
end module mod_c_nzst
