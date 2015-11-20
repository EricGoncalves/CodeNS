module mod_lp2sa
  implicit none
contains
  subroutine lp2sa( &
       v,mu,mut,dist, &
       nxn,nyn,nzn, &
       ncin,ncbd,l, &
       vol,sn,ncyc, &
       mnpar,fgam,  &
       tp,utau, &
       tn1,tn2,tn3, &
       pression,ztemp)
!
!***********************************************************************
!
!_DA  DATE_C : mai 2011 - AUTEUR:  Eric Goncalves / LEGI
!
!     ACT
!_A    Lois de paroi : approche de Smith . Modele de Spalart
!_A    - mu_t imposee a la paroi
!_A    - parois adiabatiques, 2D et 3D
!
!     INP
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
    use mod_lp2sa1
    implicit none
    integer          ::           l,       ldom,         mf,        mfb,mnpar(ip12),mnpar2(ip12)
    integer          ::  ncbd(ip41), ncin(ip41),       ncyc,         no
    double precision ::     dist(ip12),    fgam(ip42),      mu(ip12),     mut(ip12),     nxn(ip42)
    double precision ::      nyn(ip42),     nzn(ip42),pression(ip11), sn(ip31*ndir),     tn1(ip00)
    double precision ::      tn2(ip00),     tn3(ip00),      tp(ip40),    utau(ip42),  v(ip11,ip60)
    double precision ::      vol(ip11),   ztemp(ip11)
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
       if(cl(mfb)(1:3).eq.'lp4') then
!       parois adiabatiques
          call lp2sa1( &
               v,mu,mut,dist, &
               nxn,nyn,nzn, &
               ncin,ncbd,mfb,l, &
               vol,sn,ncyc,  &
               mnpar,fgam,  &
               utau, &
               tn1,tn2,tn3, &
               pression,ztemp)
!
       elseif(cl(mfb)(1:3).eq.'lp5') then
!          call lp2sa3d( &
!                v,mu,mut,dist, &
!                nxn,nyn,nzn, &
!                ncin,ncbd,mfb,l, &
!                vol,sn,ncyc,  &
!                mnpar,fgam,  &
!                utau, &
!                tn1,tn2,tn3, &
!                pression,ztemp)
       endif
    enddo
!
    return
  end subroutine lp2sa

end module mod_lp2sa
