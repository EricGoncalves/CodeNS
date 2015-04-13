module mod_lpke
  implicit none
contains
  subroutine lpke( &
       v,mu,mut,dist, &
       nxn,nyn,nzn, &
       ncin,ncbd,l, &
       mnpar,fgam,ncyc, &
       tprod,tp, &
       ztemp)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 1999-- AUTEUR : Eric Goncalves
!
!     ACT
!_A    Lois de paroi : modele k-epsilon de Jones Launder
!_A    - production de k imposee a la paroi
!_A    - epsilon imposee a la paroi
!_A    2 cas : parois adiabatiques et parois isothermes
!
!     INP
!_I    mmb        : com int (mtt       ) ; nombre de facettes d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use mod_lpke1
    use mod_lpke2
    implicit none
    integer          ::           l,       ldom,         mf,        mfb,mnpar(ip12)
    integer          ::  ncbd(ip41), ncin(ip41),       ncyc,         no
    double precision ::   dist(ip12),  fgam(ip42),    mu(ip12),   mut(ip12),   nxn(ip42)
    double precision ::    nyn(ip42),   nzn(ip42),    tp(ip40), tprod(ip00),v(ip11,ip60)
    double precision ::  ztemp(ip11)
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!
!
    nbd=0
    do no=1,mtbx
!     boucle sur toutes les frontieres
       mfb=nba(no)
       ldom=ndlb(mfb)
       if((cl(mfb)(1:2).eq.'lp').and.(l.eq.ldom)) then
!       la frontiere est une paroi et appartient au domaine en cours de traitement
          nbd=nbd+1
          lbd(nbd)=mfb
       endif
    enddo
!
    do mf=1,nbd
!     boucle sur les frontieres a traiter (parois)
       mfb=lbd(mf)
!
       if(cl(mfb)(1:3).eq.'lp2') then
!       parois adiabatiques
          call lpke1( &
               v,mu,mut,dist, &
               nxn,nyn,nzn, &
               ncin,ncbd,mfb,l, &
               mnpar,fgam,ncyc, &
               tprod,ztemp)
!
       else if(cl(mfb)(1:3).eq.'lp3') then
!       parois isothermes
          call lpke2( &
               v,mu,mut,dist, &
               nxn,nyn,nzn, &
               ncin,ncbd,mfb,l, &
               mnpar,fgam,ncyc, &
               tprod,tp,ztemp)
       endif
!
    enddo
!
!$OMP END MASTER
    return
  end subroutine lpke

end module mod_lpke
