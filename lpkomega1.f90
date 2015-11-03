module mod_lpkomega1
  implicit none
contains
  subroutine lpkomega1( &
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
!_A    Lois de paroi : modele k-omega de Wilcox
!_A    - parois adiabatiques
!_A    - production de k imposee a la paroi
!_A    - omega imposee a la paroi
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
    use proprieteflu
    use definition
    use modeleturb
    implicit none
    integer          ::   iter,     l,     m,  m0ns,    mb
    integer          ::    mfb, mnpar,  mpar,    mt,   n0c
    integer          ::     nc,  ncbd,  ncin,  ncyc,nfacns
    integer          ::     ni,   nii
    double precision ::    cmu1,   cmu2,   cmu3,     co,   dist
    double precision :: echleps,   fgam,     mu,    mup,    mut
    double precision ::      n1,     n2,     n3,    nxn,    nyn
    double precision ::     nzn,  omega,    pka,    pkb,  retur
    double precision ::     rop,     sv,     t1,     t2,     t3
    double precision ::    temp,  temp1,     tn,    top, tparoi
    double precision ::   tprod,     tt,  upyp1,    uto,      v
    double precision ::     v1t,    v1x,    v1y,    v1z,     ye
    double precision ::    yp02,     yv
    logical          :: lamin
!
!-----------------------------------------------------------------------
!
    dimension mu(ip12),mut(ip12)
    dimension nxn(ip42),nyn(ip42),nzn(ip42)
    dimension ncin(ip41),ncbd(ip41)
    dimension v(ip11,ip60),dist(ip12)
    dimension mnpar(ip12),fgam(ip42)
    dimension tprod(ip00)
    dimension temp(ip11)
!
    cmu1=1./sqrt(0.09)
    cmu2=1./(0.09**0.75)
    cmu3=0.09**0.25
    sv=110.4/tnz !air

    mt=mmb(mfb)
    m0ns=mpn(mfb)
    n0c=npc(l)
!
!       boucle sur les facettes d'une frontiere paroi
    do m=1,mt
       mb=mpb(mfb)+m
       ni=ncin(mb)
       nc=ncbd(mb)
       nfacns=m0ns+m
       nii=ni-n0c
       mpar=mnpar(ni)
!         test sur transition et regime d'ecoulement
       if((fgam(mpar).lt.1.e-3).and.(ktransi.gt.0)) then
!           laminaire
          lamin=.true.
       else
!           turbulent
          lamin=.false.
       endif
!         vitesse cellule adjacente a la paroi (cellule 1)
       v1x=v(ni,2)/v(ni,1)
       v1y=v(ni,3)/v(ni,1)
       v1z=v(ni,4)/v(ni,1)
!         normale a la paroi
       n1=nxn(nfacns)
       n2=nyn(nfacns)
       n3=nzn(nfacns)
!         tangente normee a la paroi
       tn=v1x*n1+v1y*n2+v1z*n3
       t1=v1x-tn*n1
       t2=v1y-tn*n2
       t3=v1z-tn*n3
       tt=sqrt(t1**2+t2**2+t3**2)
       t1=t1/tt
       t2=t2/tt
       t3=t3/tt
!         composante tangentielle de la vitesse dans repere paroi : v1t
       v1t=v1x*t1+v1y*t2+v1z*t3
!         temperature cellule 1 : temp1
       temp1=temp(ni)
!         temperature a la paroi : tparoi
       pka=cp*(mu(ni)/pr+mut(ni)/prt)
       pkb=mu(ni)+mut(ni)
       tparoi=temp1+0.5*pkb*v1t**2/pka
!          tparoi=temp(nc)
!         viscosite moleculaire a la paroi
       mup=mu(ni)*sqrt(tparoi/temp1)*(1.+sv/temp1)/(1.+sv/tparoi)
!         masse volumique a la paroi
       rop=v(ni,1)*temp1/tparoi
!         correction de compressibilite (loi de Van Driest)
       co=sqrt(2.*pka*tparoi/pkb)
       v1t=co*asin(v1t/co)
!         contrainte de frottement a la paroi : top
       upyp1=rop*v1t*dist(ni)/mup
       yp02=yp0**2
!         loi standard
       if(upyp1.le.yp02 .or. lamin) then
!           loi lineaire
          top=mup*v1t/dist(ni)
       else
!           loi logarithmique
          top=mup*v1t/dist(ni)
          do iter=1,10
             top=rop*v1t**2/(log(dist(ni)*sqrt(rop*top)/mup)/vkar+cllog)**2
          enddo
       endif
!        if(lamin) then
!         loi lineaire
!          top=mup*v1t/dist(ni)
!        else
!c        loi de Spalding
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
!         vitesse de frottement : uto
       uto=sqrt(top/rop)
!         calcul de la production de k -> valeur moyenne sur la cellule adjacente
       ye=2.*dist(ni)
       yv=5.*mup/(rop*uto)
       tprod(nii)=top**1.5*log(ye/yv)/(ye*vkar*sqrt(rop))
!          prodk=top**1.5*log(ye/yv)/(ye*vkar*sqrt(rop))
       if(dist(ni).lt.yv) tprod(nii)=0.
!         calcul de k en cellule 1 avec hypothese de Bradshaw
!          rok=mut(ni)*cmu1*top/(mu(ni)+mut(ni))
!          v(ni,6)=max(rok,epsk)
!         omega = k^0.5/cmu*l_eps
       retur=sqrt(rop*v(ni,6))*dist(ni)/mu(ni)
       echleps=vkar*cmu3*dist(ni)*(1.-exp(-retur/(2.*vkar*cmu2)))
       omega=sqrt(v(ni,6)/v(ni,1))/echleps
       v(ni,7)=max(v(ni,1)*omega,epse)
       v(nc,7)=v(ni,7)
!       fin boucle sur facettes paroi
    enddo
!
    return
  end subroutine lpkomega1
end module mod_lpkomega1
