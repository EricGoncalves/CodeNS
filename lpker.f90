module mod_lpker
  implicit none
contains
  subroutine lpker( &
       v,mu,mut,dist, &
       nxn,nyn,nzn, &
       ncin,ncbd,l, &
       mnpar,fgam,tprod, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       ztemp)
!
!***********************************************************************
!
!_DA  DATE_C : aout 2001 -- AUTEUR : Eric Goncalves
!
!     ACT
!_A    Lois de paroi : modele k-epsilon de Jones Launder realisable
!_A    - production de k imposee a la paroi
!_A    - epsilon imposee a la paroi
!_A    2 cas : parois adiabatiques et parois isothermes
!_A
!
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
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use mod_lpker1
    implicit none
    integer          ::           l,       ldom,         mf,        mfb,mnpar(ip12)
    integer          ::  ncbd(ip41), ncin(ip41),         no
    double precision ::   dist(ip12),  dvxx(ip00),  dvxy(ip00),  dvxz(ip00),  dvyx(ip00)
    double precision ::   dvyy(ip00),  dvyz(ip00),  dvzx(ip00),  dvzy(ip00),  dvzz(ip00)
    double precision ::   fgam(ip42),    mu(ip12),   mut(ip12),   nxn(ip42),   nyn(ip42)
    double precision ::    nzn(ip42), tprod(ip00),v(ip11,ip60), ztemp(ip11)
!
!-----------------------------------------------------------------------
!
!
!      dimension tp(ip40)
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
!
       if(cl(mfb)(1:3).eq.'lp2') then
!c      parois adiabatiques
          call lpker1( &
               v,mu,mut,dist, &
               nxn,nyn,nzn, &
               ncin,ncbd,mfb,l, &
               mnpar,fgam,tprod, &
               dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
               ztemp)
!
       else if(cl(mfb)(1:3).eq.'lp3') then
!c      parois isothermes
!          call lpker2(
!     &           v,mu,mut,dist,
!     &           nxn,nyn,nzn,
!     &           ncin,ncbd,mfb,l,
!     &           mnpar,fgam,tprod,tp
!     &           dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz)
       endif
!
    enddo
!
    return
  end subroutine lpker

end module mod_lpker
