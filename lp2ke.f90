module mod_lp2ke
  implicit none
contains
  subroutine lp2ke( &
       v,mu,mut,dist, &
       nxn,nyn,nzn, &
       ncin,ncbd,l, &
       vol,sn,ncyc, &   
       mnpar,fgam,  &   
       tprod,tp,utau, &
       tn1,tn2,tn3, & 
       pression,ztemp)
!               phipar,topx,topz)
!
!***********************************************************************
!
!_DA  DATE_C : mars 2011 - AUTEUR:  Eric Goncalves / LEGI
!
!     ACT
!_A    Lois de paroi : approche de Smith . Modele k-e de Jones-Launder
!_A    - production de k imposee a la paroi
!_A    - epsilon imposee a la paroi
!_A    2 cas : parois adiabatiques et parois isothermes
!_A
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
    use mod_lp2ke1
    implicit none
    integer          ::     l, ldom,   mf,  mfb,mnpar
    integer          ::  ncbd, ncin, ncyc,   no
    double precision ::     dist,    fgam,      mu,     mut,     nxn
    double precision ::      nyn,     nzn,pression,      sn,     tn1
    double precision ::      tn2,     tn3,      tp,   tprod,    utau
    double precision ::        v,     vol,   ztemp
!
!-----------------------------------------------------------------------
!
!
    dimension mu(ip12),mut(ip12),mnpar(ip12),dist(ip12)
    dimension nxn(ip42),nyn(ip42),nzn(ip42),fgam(ip42),utau(ip42)
    dimension ncin(ip41),ncbd(ip41),tp(ip40)
    dimension v(ip11,ip60)
    dimension tprod(ip00)
    dimension pression(ip11),ztemp(ip11),vol(ip11)     
    dimension tn1(ip00),tn2(ip00),tn3(ip00)
    dimension sn(ip31*ndir)
!      dimension phipar(ip42),topx(ip42),topz(ip42)
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
          call lp2ke1( &
               v,mu,mut,dist, &
               nxn,nyn,nzn, &
               ncin,ncbd,mfb,l, &
               vol,sn,ncyc,  &
               mnpar,fgam,  &        
               tprod,utau, &
               tn1,tn2,tn3, & 
               pression,ztemp) 
!
       elseif(cl(mfb)(1:3).eq.'lp5') then
!          call lp2ke3d( &
!                v,mu,mut,dist, &
!                nxn,nyn,nzn, &
!                ncin,ncbd,mfb,l, &
!                vol,sn,ncyc,  &
!                mnpar,fgam,   &        
!                tprod,topx,topz)  &     
       endif
!
    enddo
!      
    return
  end subroutine lp2ke

end module mod_lp2ke
