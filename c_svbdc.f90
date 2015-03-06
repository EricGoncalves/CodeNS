      subroutine c_svbdc( &
                 mot,imot,nmot, &
                 mfbi, &
                 mnc)
!
!***********************************************************************
!
!     ACT
!_A    Sauvegarde des donnees par point concernant les frontieres coincidentes
!_A    soit, pour chaque point l'indice de la maille coincidente.
!
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule coincidente
!_I    kfc        : com int              ; unite logiq, tableaux front coinc
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    mpc        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux front coinc
!
!     COM
!_C    Les donnees suvegardees ici sont celles qui peuvent etre demandees
!_C    pour l' initialisation de donnees d'une frontiere par inbdc.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
	  use boundary
	  use sortiefichier
!
!-----------------------------------------------------------------------
!
      character *32 mot(nmx)
      dimension imot(nmx)
      dimension mnc(ip43)
!
!     mfbi=nfei(mfbe)
      m0c =mpc(mfbi)
      mt  =mmb(mfbi)
!
!            call writfc( &
!                 kfc,mnc, &
!                 mt,m0c)
!
      return
      end
