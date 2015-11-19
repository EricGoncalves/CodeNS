module mod_lpke1
  implicit none
contains
  subroutine lpke1( &
       v,mu,mut,dist, &
       nxn,nyn,nzn, &
       ncin,ncbd,mfb,l, &
       mnpar,fgam,ncyc, &
       tprod,temp)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 1999-- AUTEUR : Eric Goncalves
!
!     ACT
!_A    Lois de paroi : modele k-epsilon de Jones Launder
!_A    - parois adiabatiques
!_A    - production de k imposee a la paroi
!_A    - epsilon imposee a la paroi
!_A
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
    integer          ::    mfb, mnpar(ip12),  mpar,    mt,   n0c
    integer          ::     nc,  ncbd(ip41),  ncin(ip41),  ncyc,nfacns
    integer          ::     ni,   nii
    double precision ::    cmu2, co, dist(ip12),echleps, eps
    double precision ::    fgam(ip42), mu(ip12),  mup, mut(ip12), n1
    double precision ::    n2, n3, nxn(ip42), nyn(ip42), nzn(ip42)
    double precision ::     pka,    pkb,  retur,    rop,     sv
    double precision ::      t1,     t2,     t3, temp(ip11), temp1
    double precision ::      tn,    top, tparoi, tprod(ip00), tt
    double precision ::   upyp1, uto, v(ip11,ip60), v1t, v1x
    double precision ::     v1y, v1z, ye, yp02,yv,bl,usrey
    logical          :: lamin
!
!-----------------------------------------------------------------------
!
!     lois viscosite : sutherland pour gaz et loi exponentielle pour liquide
!     loi sutherland mu=mu0*sqrt(T/T0)*(1+S/T0)/(1+S/T)
!     pour vapeur d'eau S=548K, mu0=9.73e-6 Pa.s et T0=293K
!     loi exponentielle mu=A*exp(B/T)
!     pour l'eau A=1.24e-6 Pa.s et B=1968K
      if(iflu.eq.1) then  !air
       sv=110.4/tnz
      elseif(iflu.eq.2) then !eau froide
       sv=548./tnz
       bl=1968./tnz
      elseif(iflu.eq.2) then !freon R114
!     pour vapeur R114 : S=260K, mu0=10.527e-6 Pa.s et T0=293K
!     pour le R114: A=10.336e-6 Pa.s et B=976.9738K
       sv=260./tnz
       bl=976.9738/tnz
      endif
      usrey=1./reynz
!
!   cmu1=1./sqrt(cmu)
    cmu2=1./(cmu**0.75)

    mt=mmb(mfb)
    m0ns=mpn(mfb)
    n0c=npc(l)
!
!   boucle sur les facettes d'une frontiere paroi
    do m=1,mt
       mb=mpb(mfb)+m
       ni=ncin(mb)
       nc=ncbd(mb)
       nfacns=m0ns+m
       nii=ni-n0c
       mpar=mnpar(ni)
!      test sur transition et regime d'ecoulement
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
!      normale a la paroi
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
!      tparoi=temp(nc)
!      viscosite moleculaire a la paroi
       if(iflu.eq.1) then
        mup=mu(ni)*sqrt(tparoi/temp1)*(1.+sv/temp1)/(1.+sv/tparoi)
       else
!        mup=usrey*exp(bl*(1./tparoi-1.))
        mup=mu(ni)
       endif
!      masse volumique a la paroi
       rop=v(ni,1)*temp1/tparoi
!      correction de compressibilite (loi de Van Driest)
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
!c        loi lineaire
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
!c      loi de Reichardt
!      top=mup*v1t/dist(ni)
!      do jj=1,10
!        yy=log(dist(ni)*sqrt(rop*top)/mup)
!        if (yy.lt.0.) then
!          yy=0.
!        end if
!          top=max(1.e-10,rop*v1t** 2/(2.5*log(1.+vkar*yy)+7.8*
!     &              (1.+exp(-yy/11.)-(yy/11.)*exp(-0.33*yy)))**2)
!        end do
!
!       vitesse de frottement : uto
       uto=sqrt(top/rop)
!       calcul de la production de k -> valeur moyenne sur la cellule adjacente
       ye=2.*dist(ni)
       yv=5.*mup/(rop*uto)
       tprod(nii)=top**1.5*log(ye/yv)/(ye*vkar*sqrt(rop))
!     &         + uto*dpdx(nii)*(ye-yv)/(vkar*ye) &
!               + top*dvdx*(ye-yv)/ye &
!               + dpdx(nii)*dvdx*0.5*(ye**2-yv**2)/ye
       if(dist(ni).lt.yv) tprod(nii)=0.
!       k et epsilon de Mohammadi
!       rok=top*cmu1*(min(yplus1/100.,1.))**2
!       v(ni,6)=max(rok,epsk)
!       epsilon = k^1.5/l_eps
       retur=sqrt(v(ni,1)*v(ni,6))*dist(ni)/mu(ni)
       echleps=vkar*cmu2*dist(ni)*(1.-exp(-retur/(2.*vkar*cmu2)))
       eps=(v(ni,6)/v(ni,1))**1.5/echleps
       v(ni,7)=max(v(ni,1)*eps,epse)
       v(nc,7)=0.
!     fin boucle sur facettes paroi
    enddo
!
    return
  end subroutine lpke1
end module mod_lpke1
