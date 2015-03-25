module mod_lpsa2
  implicit none
contains
  subroutine lpsa2( &
       v,mu,mut,dist, &
       nxn,nyn,nzn, &
       ncin,ncbd,mfb,l, &
       mnpar,fgam,ncyc,tp, &
       temp)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 1999-- AUTEUR :  Eric Goncalves
!
!     ACT
!_A    Lois de paroi : modele de Spalart&Allmaras
!_A    - parois isothermes
!_A    - nut imposee a la paroi
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
  double precision ::           c1,          c2,          c3,          c4,          ca
  double precision ::           cb,         cta,         ctb,       denom,  dist(ip12)
  double precision ::         dudy,  fgam(ip42),    mu(ip12),         mup,   mut(ip12)
  double precision ::           n1,          n2,          n3,   nxn(ip42),   nyn(ip42)
  double precision ::    nzn(ip42),         rc4,    rnutilde,         rop,          sv
  double precision ::           t1,          t2,          t3,  temp(ip11),       temp1
  double precision ::           tn,         top,    tp(ip40),          tt,       upyp1
  double precision :: v(ip11,ip60),         v1t,         v1x,         v1y,         v1z
  double precision ::         yp02,      yplus1
  logical          :: lamin
!
!-----------------------------------------------------------------------
!
!
!
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
       mpar=mnpar(ni)
       nii=ni-n0c
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
!         masse volumique a la paroi
       rop=v(ni,1)*temp1/tp(m)
!         viscosite moleculaire a la paroi
       mup=mu(ni)*sqrt(tp(m)/temp1)*(1.+sv/temp1)/(1.+sv/tp(m))
!         correction de compressibilite (loi de Van Driest)
       cta=(mu(ni)+mut(ni))/(cp*(mu(ni)/pr+mut(ni)/prt))
       ctb=cta/(2.*tp(m))
       denom = (tp(m)-temp1-cta*0.5*v1t**2)*sqrt(temp1/tp(m)) + &
            tp(m)-temp1+cta*0.5*v1t**2
       v1t=(1./sqrt(ctb))*asin(2.*sqrt(ctb)*(tp(m)-temp1)*v1t/denom)
!         contrainte de frottement a la paroi : top
       upyp1=rop*v1t*dist(ni)/mup
       yp02=yp0**2
!         loi standard
       if(upyp1.le.yp02 .or. lamin) then
!           loi lineaire
          top=mup*v1t/dist(ni)
          dudy=top/mup
       else
!           loi logarithmique
          top=mup*v1t/dist(ni)
          do iter=1,10
             top=rop*v1t**2/(log(dist(ni)*sqrt(rop*top)/mup)/vkar+cllog)**2
          enddo
          dudy=sqrt(top/rop)/(vkar*dist(ni))
       endif
!
!        loi raffinee avec interpolation
!         if(upyp1.le.9.) then
!          loi lineaire
!           top=max(1.e-10,mup*v1t/dist(ni))
!       else
!         yplus0=3.
!         do ii=1,10
!           yplus0=upyp1/(log(yplus0)/vkar+cllog)
!          end do
!           if(yplus0.gt.40.) then
!            loi logarithmique
!           top=max(1.e-10,mup*v1t/dist(ni))
!           do iter=1,10
!             top=max(1.e-10,rop)*v1t**2/(log(dist(ni)*sqrt(rop
!     &               *top)/mup)/vkar+cllog)**2)
!            end do
!           else
!            region tampon : interpolation par polynome de degre 4
!           top=max(1.e-10,mup*v1t/dist(ni))
!           do jj=1,10
!             yy=log(dist(ni)*sqrt(rop*top)/mup)
!              top=max(1.e-10,rop*v1t** 2/(0.17962*yy**4-
!     &                  2.2117*yy**3+9.2052*yy**2-10.804*yy
!     &             +6.4424)**2)
!            end do
!           end if
!        end if
!
!         yplus cellule 1
       yplus1=dist(ni)*sqrt(top*rop)/mup
!         calcul de nu_tilde en cellule 1
       mut(ni)=v(ni,1)*vkar**2*dist(ni)**2*dudy* &
            (1.-exp(-yplus1/26.))**2
       ca=mut(ni)
       cb=mu(ni)**3*mut(ni)*cv1**3
       c1=cb*sqrt(3.)*sqrt(256.*cb+27.*ca**4)
       c2=0.5*cb*ca**2+(1./18.)*c1
       c3=-0.5*cb*ca**2+(1./18.)*c1
       c4=4.*abs(c3**(1./3.)-c2**(1./3.))
       rc4=sqrt(ca**2+c4)
       rnutilde=0.25*(ca+rc4)+0.25*sqrt(2.)*sqrt(ca**2*rc4-2.* &
            c3**(1./3.)*rc4+2.*c2**(1./3.)*rc4+ca**3)/rc4**0.5
       v(ni,6)=max(rnutilde,epsk)
       v(nc,6)=0.
!         ATTENTION! si probleme en maillage fin, utiliser nu=kappa*dist*utau
!          rnutilde=v(ni,1)*kappa*dist(ni)*sqrt(top/rop)
!
!       fin boucle sur facettes paroi
    end do
!
    return
  end subroutine lpsa2

end module mod_lpsa2
