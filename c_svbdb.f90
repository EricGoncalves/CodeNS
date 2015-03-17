module mod_c_svbdb
  implicit none
contains
  subroutine c_svbdb( &
       mot,imot,nmot, &
       mfbi, &
       ncin)
!
!***********************************************************************
!
!     ACT
!_A    Sauvegarde des donnees de base par point d'une frontiere
!_A    soit, l'indice du centre de la maille qui s'appuie sur la
!_A    frontiere.
!
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    kfb        : com int              ; unite logiq, tableaux de base front
!
!     COM
!_C    Les donnees suvegardees ici sont celles qui peuvent etre demandees
!_C    pour l' initialisation des donnees de base d'une frontiere par inbdb.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use boundary
    implicit none
    integer          :: imot,  m0,mfbi,  mt,ncin
    integer          :: nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    dimension imot(nmx)
    dimension ncin(ip41)
!
!     mfbi=nfei(mfbe)
    m0  =mpb(mfbi)
    mt  =mmb(mfbi)
!
!      call writfb( &
!                 kfb,ncin, &
!                 mt,m0)
!
    return
  end subroutine c_svbdb
end module mod_c_svbdb
