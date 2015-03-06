      subroutine c_svbd( &
                 mot,imot,nmot, &
                 mfbi, &
                 ncbd)
!
!***********************************************************************
!
!     ACT
!_A    Sauvegarde des donnees necessaires a l'identification d'une frontiere
!_A    soit, son nombre de points et le reperage de ces points dans les
!_A    tableaux domaines.
!
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    kfi        : com int              ; unite logiq, indices pts frontieres
!
!     COM
!_C    Les donnees suvegardees ici sont celles qui peuvent etre demandees pour
!_C    la creation d'une frontiere par crbd.
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
      dimension ncbd(ip41)
!
      m0  =mpb(mfbi)
      mt  =mmb(mfbi)
!
!       call writfi( &
!                 kfi,ncbd, &
!                 mt,m0)
!
      return
      end
