module mod_lpke2
  implicit none
contains
  subroutine lpke2( &
       v,mu,mut,dist, &
       nxn,nyn,nzn, &
       ncin,ncbd,mfb,l, &
       mnpar,fgam,ncyc, &
       tprod,tp,temp)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 1999-- AUTEUR :  Eric Goncalves
!
!     ACT
!_A    Lois de paroi : modele k-epsilon de Jones Launder
!_A    - parois isothermes
!_A    - production de k imposee a la paroi
!_A    - epsilon imposee a la paroi
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
    double precision ::         cmu2,         cta,         ctb,         ctc,         cts
    double precision ::   dist(ip12),     echleps,         eps,  fgam(ip42),    mu(ip12)
    double precision ::          mup,   mut(ip12),          n1,          n2,          n3
    double precision ::    nxn(ip42),   nyn(ip42),   nzn(ip42),       prodk,       retur
    double precision ::           sv,          t1,          t2,          t3,  temp(ip11)
    double precision ::        temp1,          tn,         top,    tp(ip40), tprod(ip00)
    double precision ::           tt,       upyp1,         uto,v(ip11,ip60),         v1t
    double precision ::          v1x,         v1y,         v1z,          ye,        yp02
    double precision ::           yv
    logical          :: lamin
!
!-----------------------------------------------------------------------
!
!
!
!     ctl=vkar/(cmu**0.75)
!      cmu1=1./sqrt(cmu)
    cmu2=1./(cmu**0.75)
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
       if((fgam(mpar).lt.1.e-3).and.(ktransi.gt.0)) then
!         laminaire
          lamin=.true.
       else
!         turbulent
          lamin=.false.
       end if
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
!       masse volumique a la paroi
       v(nc,1)=v(ni,1)*temp1/tp(m)
!       viscosite moleculaire a la paroi
       mup=mu(ni)*sqrt(tp(m)/temp1)*(1.+sv/temp1)/(1.+sv/tp(m))
!       correction de compressibilite (loi de Van Driest)
       cta=(mu(ni)+mut(ni))/(cp*(mu(ni)/pr+mut(ni)/prt))
       ctb=cta/(2.*tp(m))
       ctc=(tp(m)-temp1-0.5*cta*v1t**2)/(tp(m)*v1t)
       cts=sqrt(4.*ctb+ctc**2)
       v1t=(1./sqrt(ctb))*(asin((2.*ctb*v1t+ctc)/cts)-asin(ctc/cts))
!       contrainte de frottement a la paroi : top
       upyp1=v(nc,1)*v1t*dist(ni)/mup
       yp02=yp0**2
!        if(upyp1.le.yp02) then
       if(upyp1.le.yp02 .or. lamin) then
!         loi lineaire
          top=mup*v1t/dist(ni)
       else
!         loi logarithmique
          top=mup*v1t/dist(ni)
          do iter=1,10
             top=v(nc,1)*v1t**2/(log(dist(ni)* &
                  sqrt(v(nc,1)*top)/mup)/vkar+cllog)**2
          end do
       end if
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
!c           loi logarithmique
!           top=max(1.e-10,mup*v1t/dist(ni))
!           do iter=1,10
!             top=max(1.e-10,v(nc,1)*v1t**2/(log(dist(ni)*sqrt(v(nc,1)
!     &               *top)/mup)/vkar+cllog)**2)
!            end do
!           else
!c           region tampon : interpolation par polynome de degre 4
!           top=max(1.e-10,mup*v1t/dist(ni))
!           do jj=1,10
!             yy=log(dist(ni)*sqrt(v(nc,1)*top)/mup)
!              top=max(1.e-10,v(nc,1)*v1t** 2/(0.17962*yy**4-
!     &                  2.2117*yy**3+9.2052*yy**2-10.804*yy
!     &             +6.4424)**2)
!            end do
!           end if
!       end if
!      vitesse de frottement : uto
       uto=sqrt(top/v(nc,1))
!       calcul de la production de k -> valeur moyenne sur la cellule adjacente
       ye=2.*dist(ni)
       yv=5.*mup/(v(nc,1)*uto)
       prodk=top**1.5*log(ye/yv)/(ye*vkar*sqrt(v(nc,1)))
!     &           + top*dvdx*(ye-yv)/ye
!       calcul de k en cellule 1 avec hypothese de Bradshaw
!       rok=mut(ni)*cmu1*top/(mu(ni)+mut(ni))
!        v(ni,6)=max(rok,epsk)
!         k et epsilon de Mohammadi
!         epsilon = k^1.5/l_eps
       retur=sqrt(v(ni,1)*v(ni,6))*dist(ni)/mu(ni)
       echleps=vkar*cmu2*dist(ni)*(1.-exp(-retur/(2.*vkar*cmu2)))
       eps=(v(ni,6)/v(ni,1))**1.5/echleps
       v(ni,7)=max(v(ni,1)*eps,epse)
       v(nc,7)=0.
       tprod(nii)=prodk
       if(dist(ni).lt.yv) tprod(nii)=0.
!     fin boucle sur facettes paroi
    enddo
!
    return
  end subroutine lpke2
end module mod_lpke2
