module mod_lpkl1
  implicit none
contains
  subroutine lpkl1( &
       v,mu,mut,dist, &
       nxn,nyn,nzn, &
       ncin,ncbd,mfb,l, &
       mnpar,fgam, &
       tprod,ncyc, &
       temp)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 1999-- AUTEUR :  Eric Goncalves
!
!     ACT
!_A    Lois de paroi, modele k-l de Smith, parois adiabatiques
!_A    - la production de k est imposee a la paroi
!_A    - la longuer l est imposee a la paroi
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
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use proprieteflu
    use definition
    use modeleturb
    implicit none
    integer          ::        iter,          l,          m,       m0ns,         mb
    integer          ::         mfb,mnpar(ip12),       mpar,         mt,        n0c
    integer          ::          nc, ncbd(ip41), ncin(ip41),       ncyc,     nfacns
    integer          ::          ni,        nii
    double precision ::         cmu1,          co,  dist(ip12),  fgam(ip42),    mu(ip12)
    double precision ::          mup,   mut(ip12),          n1,          n2,          n3
    double precision ::    nxn(ip42),   nyn(ip42),   nzn(ip42),         pka,         pkb
    double precision ::         rhol,         rop,          sv,          t1,          t2
    double precision ::           t3,  temp(ip11),       temp1,          tn,         top
    double precision ::       tparoi, tprod(ip00),          tt,       upyp1,         uto
    double precision :: v(ip11,ip60),         v1t,         v1x,         v1y,         v1z
    double precision ::           ye,        yp02,          yv
    logical          :: lamin
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!
!
    cmu1=1./sqrt(cmu)
    sv=110.4/tnz !air

    mt=mmb(mfb)
    m0ns=mpn(mfb)
    n0c=npc(l)
!
!     boucle sur les facettes d'une frontiere paroi
    do m=1,mt
       mb=mpb(mfb)+m
       ni=ncin(mb)
       nc=ncbd(mb)
       nfacns=m0ns+m
       nii=ni-n0c
       mpar=mnpar(ni)
!       test sur transition et regime d'ecoulement
       if((fgam(mpar).lt.1.e-3).and.(ktransi.gt.0)) then
!         laminaire
          lamin=.true.
       else
!         turbulent
          lamin=.false.
       endif
!       vitesse cellule adjacente a la paroi (cellule 1)
       v1x=v(ni,2)/v(ni,1)
       v1y=v(ni,3)/v(ni,1)
       v1z=v(ni,4)/v(ni,1)
!       normale a la paroi
       n1=nxn(nfacns)
       n2=nyn(nfacns)
       n3=nzn(nfacns)
!       tangente normee a la paroi
       tn=v1x*n1+v1y*n2+v1z*n3
       t1=v1x-tn*n1
       t2=v1y-tn*n2
       t3=v1z-tn*n3
       tt=sqrt(t1**2+t2**2+t3**2)
       t1=t1/tt
       t2=t2/tt
       t3=t3/tt
!       composante tangentielle de la vitesse dans repere paroi : v1t
       v1t=v1x*t1+v1y*t2+v1z*t3
!       temperature cellule 1 : temp1
       temp1=temp(ni)
!       temperature a la paroi : tparoi
       pka=cp*(mu(ni)/pr+mut(ni)/prt)
       pkb=mu(ni)+mut(ni)
       tparoi=temp1+0.5*pkb*v1t**2/pka
!        tparoi=temp(nc)
!       viscosite moleculaire a la paroi
       mup=mu(ni)*sqrt(tparoi/temp1)*(1.+sv/temp1)/(1.+sv/tparoi)
!       masse volumique a la paroi
       rop=v(ni,1)*temp1/tparoi
!       correction de compressibilite (loi de Van Driest)
       co=sqrt(2.*pka*tparoi/pkb)
       v1t=co*asin(v1t/co)
!       contrainte de frottement a la paroi : top
       upyp1=rop*v1t*dist(ni)/mup
       yp02=yp0**2
       if(upyp1.le.yp02 .or. lamin) then
!         loi lineaire
          top=mup*v1t/dist(ni)
       else
!         loi logarithmique
          top=mup*v1t/dist(ni)
          do iter=1,10
             top=rop*v1t**2/(log(dist(ni)*sqrt(rop*top)/mup)/vkar+cllog)**2
          enddo
       endif
!
!        if(lamin) then
!         loi lineaire
!          top=mup*v1t/dist(ni)
!        else
!         loi de Spalding
!          top=mup*v1t/dist(ni)
!          upl=v1t*sqrt(rop/top)
!          do jj=1,15
!            fu=rop*dist(ni)*v1t/(mup*upl) - upl - exp(-vkar*cllog)*
!     &         (exp(vkar*upl) -1. -vkar*upl - 0.5*(vkar*upl)**2 -
!     &          ((vkar*upl)**3)/6.)
!            dfu=-rop*dist(ni)*v1t/(mup*upl**2)-1.-exp(-vkar*cllog)*
!     &          vkar*(exp(vkar*upl) -1. -vkar*upl - 0.5*(vkar*upl)**2)
!            upl=upl-fu/dfu
!          enddo
!        endif
!        top=rop*(v1t/upl)**2
!
!       loi de Reichardt
!        top=max(1.e-10,mup*v1t/dist(ni))
!        do jj=1,10
!         yy=log(dist(ni)*sqrt(rop*top)/mup)
!          top=max(1.e-10,rop*v1t** 2/(2.5*log(1.+vkar*yy)+7.8*
!     &        (1.+exp(-yy/11.)-(yy/11.)*exp(-0.33*yy)))**2)
!        end do
!
!       calcul de la production de k -> valeur moyenne sur la cellule adjacente
       uto=sqrt(top/rop)
       ye=2.*dist(ni)
       yv=5.*mup/(rop*uto)
!         yplus1=rop*uto*dist(ni)/mup
       tprod(nii)=top**1.5*log(ye/yv)/(ye*vkar*sqrt(rop))
!     &             + top*dvdx*(ye-yv)/ye
       if(dist(ni).lt.yv) tprod(nii)=0.
!       calcul de k en cellule 1 avec hypothese de Bradshaw
!        rok=mut(ni)*cmu1*top/(mu(ni)+mut(ni))
!        v(ni,6)=max(rok,epsk)
!       calcul de l en cellule 1
       rhol=vkar*dist(ni)*v(ni,1)
       v(ni,7)=max(rhol,epse)
       v(nc,7)=0.
!       fin boucle sur facettes paroi
    enddo
!
!$OMP END MASTER
    return
  end subroutine lpkl1

end module mod_lpkl1
