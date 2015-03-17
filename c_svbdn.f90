module mod_c_svbdn
  implicit none
contains
  subroutine c_svbdn( &
       mot,imot,nmot, &
       mfbi, &
       nxn,nyn,nzn)
!
!***********************************************************************
!
!     ACT
!_A    Sauvegarde des donnees par point d'une frontiere a normales
!_A    stockees soit, les trois composantes de la normale en chaque
!_A    point.
!
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    kfn        : com int              ; unite logiq, tableaux normales
!
!     COM
!_C    Les donnees suvegardees ici sont celles qui peuvent etre demandees
!_C    pour l' initialisation de donnees d'une frontiere par inbdn.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use boundary
    implicit none
    integer          :: imot, m0n,mfbi,  mt,nmot
    double precision :: nxn,nyn,nzn
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    dimension nxn(ip42),nyn(ip42),nzn(ip42)
    dimension imot(nmx)
!
    m0n =mpn(mfbi)
    mt  =mmb(mfbi)
!
!            call writfn( &
!                 kfn,nxn,nyn,nzn, &
!                 mt,m0n)
!
    return
  end subroutine c_svbdn
end module mod_c_svbdn
