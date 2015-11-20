module mod_lpkomega
  implicit none
contains
  subroutine lpkomega( &
       v,mu,mut,dist, &
       nxn,nyn,nzn, &
       ncin,ncbd,l, &
       mnpar,fgam, &
       tprod,ncyc,tp, &
       ztemp)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 1999-- AUTEUR : Eric Goncalves
!
!     ACT
!_A    Lois de paroi : modele k-omega de Wilcox
!_A    - production de k imposee a la paroi
!_A    - omega imposee a la paroi
!_A    2 cas : parois adiabatiques et parois isothermes
!_A
!
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    u          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    cl         : com char(mtb       ) ; type de cond lim a appliquer
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    mmb        : com int (mtt       ) ; nombre de facettes d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!_I    mnpar      : arg real(ip12      ) ; pointeur dans tableaux front normales
!_I                                        stockees du point de rattach normale
!_I    fgam       : arg real(ip42      ) ; fonction d'intermittence pour
!_I                                        transition
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use mod_lpkomega1
    use mod_lpkomega2
    implicit none
    integer          ::           l,       ldom,         mf,        mfb,mnpar(ip12),mnpar2(ip12)
    integer          ::  ncbd(ip41), ncin(ip41),       ncyc,         no
    double precision ::   dist(ip12),  fgam(ip42),    mu(ip12),   mut(ip12),   nxn(ip42)
    double precision ::    nyn(ip42),   nzn(ip42),    tp(ip40), tprod(ip00),v(ip11,ip60)
    double precision ::  ztemp(ip11)
!
!-----------------------------------------------------------------------
!
!
!
    nbd=0
    do no=1,mtbx
!c    boucle sur toutes les frontieres
       mfb=nba(no)
       ldom=ndlb(mfb)
       if((cl(mfb)(1:2).eq.'lp').and.(l.eq.ldom)) then
!c      la frontiere est une paroi et appartient au domaine en cours de traitement
          nbd=nbd+1
          lbd(nbd)=mfb
       endif
    enddo
!
    do mf=1,nbd
!c    boucle sur les frontieres a traiter (parois)
       mfb=lbd(mf)
       if(cl(mfb)(1:3).eq.'lp2') then
!c      parois adiabatiques
          call lpkomega1( &
               v,mu,mut,dist, &
               nxn,nyn,nzn, &
               ncin,ncbd,mfb,l, &
               mnpar,fgam, &
               tprod,ncyc, &
               ztemp)
!
       else if(cl(mfb)(1:3).eq.'lp3') then
!c      parois isothermes
          call lpkomega2( &
               v,mu,mut,dist, &
               nxn,nyn,nzn, &
               ncin,ncbd,mfb,l, &
               mnpar,fgam, &
               tprod,ncyc,tp, &
               ztemp)
       endif
!
    enddo
!
    return
  end subroutine lpkomega
end module mod_lpkomega
