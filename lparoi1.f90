module mod_lparoi1
  implicit none
contains
  subroutine lparoi1( &
       lm,ncyc, &
       v,mu,mut,dist, &
       nxn,nyn,nzn, &
       ncin,ncbd,mfbm, &
       toxx,toxy,toxz,toyy,toyz,tozz, &
       qcx,qcy,qcz, &
       mnpar,fgam, &
       temp)
!
!***********************************************************************
!
!_DA  DATE_C : decembre 1998-- AUTEUR : Eric GONCALVES
!
!     ACT
!_A    Lois de paroi analytiques
!_A    - parois adiabatiques
!_A    - tenseur des contraintes imposee a la paroi
!_A    - vecteur flux de chaleur imposee a la paroi
!_A
!_A    Les densites de flux numeriques a la paroi sont evaluees a partir
!_A    du champ aerodynamique au niveau des cellules adjacentes aux parois.
!
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    npfb       : com int (lt        ) ; pointeur fin de domaine precedent
!_I                                        dans tableau toutes facettes
!_I    nnn        : com int (lt        ) ; nombre de noeuds du domaine (dont fictif)
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_I                                        contraintes visqueuses
!_I    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_I                                        contraintes visqueuses
!_I    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_I                                        contraintes visqueuses
!_I    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_I                                        contraintes visqueuses
!_I    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_I                                        contraintes visqueuses
!_I    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_I                                        contraintes visqueuses
!_I    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_I    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_I    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
!_I    cl         : com char(mtb       ) ; type de cond lim a appliquer
!_I    img        : com int              ; niveau de grille (multigrille)
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    mmb        : com int (mtt       ) ; nombre de facettes d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    npfb       : com int (lt        ) ; pointeur fin de domaine precedent
!_I                                        dans tableau toutes facettes
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
    use modeleturb
    use definition
    implicit none
    integer          ::        iter,         lm,          m,       m0ns,         mb
    integer          ::        mfbm,mnpar(ip12),       mpar,         mt,        n0c
    integer          ::          nc, ncbd(ip41), ncin(ip41),       ncyc,     nfacns
    integer          ::          ni,        nii
    double precision ::           co,  dist(ip12),  fgam(ip42),    mu(ip12),         mup
    double precision ::    mut(ip12),          n1,          n2,          n3,   nxn(ip42)
    double precision ::    nyn(ip42),   nzn(ip42),         pka,         pkb,         qc1
    double precision ::    qcx(ip12),   qcy(ip12),   qcz(ip12),         rop,          sv
    double precision ::           t1,          t2,          t3,  temp(ip11),       temp1
    double precision ::           tn,         top,  toxx(ip12),  toxy(ip12),  toxz(ip12)
    double precision ::   toyy(ip12),  toyz(ip12),  tozz(ip12),      tparoi,          tt
    double precision ::        upyp1,v(ip11,ip60),         v1t,         v1x,         v1y
    double precision ::          v1z,        yp02
    logical          :: lamin
!
!-----------------------------------------------------------------------
!
!
!     lois viscosite : sutherland pour gaz et loi exponentielle pour liquide
!     loi sutherland mu=mu0*sqrt(T/T0)*(1+S/T0)/(1+S/T)
!     pour vapeur d'eau S=548K, mu0=9.73e-6 Pa.s et T0=293K
!     loi exponentielle mu=A*exp(B/T)
!     pour l'eau A=1.24e-6 Pa.s et B=1968K
!      if(iflu.eq.1) then  !eau
!        sv=548./tnz
!        s0=tnz/293.
!        mu0=9.73E-6
!        al=1.214E-6
!        bl=1968./tnz
!      elseif(iflu.eq.2) then !freon R114
!     pour vapeur R114 : S=260K, mu0=10.527e-6 Pa.s et T0=293K
!     pour le R114: A=10.336e-6 Pa.s et B=976.9738K
!        sv=260./tnz
!        s0=tnz/293.
!        mu0=10.527E-6
!        al=10.336E-6
!        bl=976.9738/tnz
!      elseif(iflu.eq.3) then !air
!       sv=110.4/tnz
!        STOP
!      endif
!
    sv=110.4/tnz !air
    mt=mmb(mfbm)
    m0ns=mpn(mfbm)
    n0c=npc(lm)
!
!     boucle sur les facettes d'une frontiere paroi
    do m=1,mt
       mb=mpb(mfbm)+m
       ni=ncin(mb)
       nc=ncbd(mb)
       nfacns=m0ns+m
       mpar=mnpar(ni)
       nii=ni-n0c
!       test sur transition et regime d'ecoulement
       if((fgam(mpar).lt.1.e-3).and.(ktransi.gt.0)) then
!         laminaire
          lamin=.true.
       else
!         turbulent
          lamin=.false.
       endif
!      vitesse cellule adjacente a la paroi (cellule 1)
       v1x=v(ni,2)/v(ni,1)
       v1y=v(ni,3)/v(ni,1)
       v1z=v(ni,4)/v(ni,1)
!       normale a la paroi
       n1=nxn(nfacns)
       n2=nyn(nfacns)
       n3=nzn(nfacns)
!      tangente normee a la paroi
       tn=v1x*n1+v1y*n2+v1z*n3
       t1=v1x-tn*n1
       t2=v1y-tn*n2
       t3=v1z-tn*n3
       tt=sqrt(t1**2+t2**2+t3**2)
       t1=t1/tt
       t2=t2/tt
       t3=t3/tt
!      composante tangentielle de la vitesse dans repere paroi : v1t
       v1t=v1x*t1+v1y*t2+v1z*t3
!      temperature cellule 1 : temp1
       temp1=temp(ni)
!      temperature a la paroi : tparoi
       pka=cp*(mu(ni)/pr+mut(ni)/prt)
       pkb=mu(ni)+mut(ni)
       tparoi=temp1+0.5*pkb*v1t**2/pka
!      masse volumique a la paroi
       rop=v(ni,1)*temp1/tparoi
!      viscosite moleculaire a la paroi
       mup=mu(ni)*sqrt(tparoi/temp1)*(1.+sv/temp1)/(1.+sv/tparoi)
!      correction de compressibilite (loi de Van Driest)
       co=sqrt(2.*pka*tparoi/pkb)
       v1t=co*asin(v1t/co)
!       contrainte de frottement a la paroi : top
       upyp1=rop*v1t*dist(ni)/mup
       yp02=yp0**2
       if(upyp1.le.yp02 .or. lamin) then
!        loi lineaire
          top=mup*v1t/dist(ni)
       else
!        loi logarithmique
          top=mup*v1t/dist(ni)
          do iter=1,10
             top=rop*v1t**2/(log(dist(ni)*sqrt(rop*top)/mup)/vkar+cllog)**2
          enddo
       endif
!
!        if(lamin) then
!        loi lineaire
!          top=mup*v1t/dist(ni)
!        else
!        loi de Spalding
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
!       tenseur des contraintes a la paroi dans repere general
       toxx(nc)=2*t1*n1*top
       toyy(nc)=2*t2*n2*top
       tozz(nc)=2*t3*n3*top
       toxy(nc)=(t2*n1+t1*n2)*top
       toxz(nc)=(t1*n3+t3*n1)*top
       toyz(nc)=(t2*n3+t3*n2)*top
!       tenseur des contraintes en 1 dans repere general
       toxx(ni)=toxx(nc)
       toyy(ni)=toyy(nc)
       tozz(ni)=tozz(nc)
       toxy(ni)=toxy(nc)
       toxz(ni)=toxz(nc)
       toyz(ni)=toyz(nc)
!       flux de chaleur dans cellule 1 dans repere general
!       ATTENTION! le signe 'moins' provient de la convention utilisee dans le code
       qc1=top*((v(ni,2)*t1+v(ni,3)*t2+v(ni,4)*t3)/v(ni,1))
       qcx(ni)=-n1*qc1
       qcy(ni)=-n2*qc1
       qcz(ni)=-n3*qc1
!       flux de chaleur a la paroi dans repere general
       qcx(nc)=0.
       qcy(nc)=0.
       qcz(nc)=0.
!
!     fin boucle sur facettes d'une frontiere paroi
    enddo
!
    return
  end subroutine lparoi1
end module mod_lparoi1
