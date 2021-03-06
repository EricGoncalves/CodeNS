module mod_sch_ausmp_pond
  implicit none
contains
  subroutine sch_ausmp_pond( &
       lm,ityprk, &
       u,v,d,ff, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       equat, &
       sn,lgsnlt, &
       rhol,ul,vl,wl,pl,rhor,ur,vr,wr,prr, &
       ps, &
       cvi,cvj,cvk, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DA  DATE_C : avril 2002 - Eric GONCALVES / Sinumef
!
!     ACT
!_A    Calcul des bilans de flux physiques: Schema AUSM+ de Liou.
!_A    Ordres 1, 2 et 3 avec extrapolation MUSCL.
!_A    Ponderation pour maillage irregulier.
!_A    Limiteur de pentes.
!
!***********************************************************************
!
    use para_var
    use para_fige
    use maillage
    use schemanum
    use proprieteflu
    implicit none
    integer          ::     i1,  i1m1,  i1p1,    i2
    integer          ::   i2m1,  ind1,  ind2,ityprk
    integer          ::      j,    j1,  j1m1,  j1p1,    j2
    integer          ::   j2m1,     k,    k1,  k1m1
    integer          ::   k1p1,    k2,  k2m1,  kdir
    integer          :: lgsnlt,    lm,     m,     n,   n0c
    integer          ::     n1,   nci,   ncj,   nck,   nid
    integer          ::   nijd,  ninc,   njd
    double precision ::                   ai,                  al,                  ar
    double precision ::                  cal,                 car,                 cbl,                 cbr,                 cmi
    double precision ::                  cml,                 cmr,         cmui1(ip21),         cmui2(ip21),         cmuj1(ip21)
    double precision ::          cmuj2(ip21),         cmuk1(ip21),         cmuk2(ip21),                cnds,           cvi(ip21)
    double precision ::            cvj(ip21),           cvk(ip21),        d(ip11,ip60),                  dm,                 fex
    double precision ::                  fey,                 fez,       ff(ip11,ip60),                 fxx,                 fxy
    double precision ::                  fxz,                 fyy,                 fyz,                 fzz,                 hi1
    double precision ::                  hi2,                 hi3,                 hi4,                 hi5,                 hj1
    double precision ::                  hj2,                 hj3,                 hj4,                 hj5,                 hk1
    double precision ::                  hk2,                 hk3,                 hk4,                 hk5,                  hl
    double precision ::                   hr,            pl(ip00),           prr(ip00),            ps(ip11),                 psi
    double precision ::            qcx(ip12),           qcy(ip12),           qcz(ip12),          rhol(ip00)
    double precision ::           rhor(ip00),                 si0,                 si1,                 si2,                 si3
    double precision ::                  si4,                 sj0,                 sj1,                 sj2,                 sj3
    double precision ::                  sj4,                 sk0,                 sk1,                 sk2,                 sk3
    double precision ::                  sk4,sn(lgsnlt,nind,ndir),          toxx(ip12),          toxy(ip12)
    double precision ::           toxz(ip12),          toyy(ip12),          toyz(ip12),          tozz(ip12),        u(ip11,ip60)
    double precision ::             ul(ip00),            ur(ip00),        v(ip11,ip60),            vl(ip00),            vr(ip00)
    double precision ::             wl(ip00),            wr(ip00)
    double precision,allocatable :: r1(:),r2(:),r3(:),r4(:),r5(:)
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
!
    ALLOCATE(r1(ip00),r2(ip00),r3(ip00),r4(ip00),r5(ip00))

    n0c=npc(lm)
    i1=ii1(lm)
    i2=ii2(lm)
    j1=jj1(lm)
    j2=jj2(lm)
    k1=kk1(lm)
    k2=kk2(lm)

    nid = id2(lm)-id1(lm)+1
    njd = jd2(lm)-jd1(lm)+1
    nijd = nid*njd

    i1p1=i1+1
    j1p1=j1+1
    k1p1=k1+1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1

    nci = inc(1,0,0)
    ncj = inc(0,1,0)
    nck = inc(0,0,1)
!
!-----calcul des densites de flux visqueux--------------------------------
!
    if(equat(3:5).eq.'2dk') then
       ind1 = indc(i1m1,j1m1,k1  )
       ind2 = indc(i2  ,j2  ,k2m1)
    elseif(equat(3:4).eq.'3d') then
       ind1 = indc(i1m1,j1m1,k1m1)
       ind2 = indc(i2  ,j2  ,k2  )
    endif
    do n=ind1,ind2
       m=n-n0c
       u(n,1)=0.
       u(n,2)=0.
       u(n,3)=0.
       u(n,4)=0.
       u(n,5)=0.
       d(n,1)=(toxx(n)*v(n,2)+toxy(n)*v(n,3)+toxz(n)*v(n,4))/v(n,1)+qcx(n)
       d(n,2)=(toyy(n)*v(n,3)+toxy(n)*v(n,2)+toyz(n)*v(n,4))/v(n,1)+qcy(n)
       d(n,3)=(tozz(n)*v(n,4)+toxz(n)*v(n,2)+toyz(n)*v(n,3))/v(n,1)+qcz(n)
    enddo
!
!*********************************************************************
!      calcul des flux numeriques par direction
!*********************************************************************
!
!------direction i----------------------------------------------
!
    kdir=1
    ninc=nci
!
!-----definition des variables extrapolees------------------------
!
    if(ilim.eq.1) then
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1,j,k)
             ind2 = indc(i2,j,k)
             do n=ind1,ind2
                m=n-n0c
                r1(m)=v(n,1)-v(n-ninc,1)
                r2(m)=v(n,2)/v(n,1)-v(n-ninc,2)/v(n-ninc,1)
                r3(m)=v(n,3)/v(n,1)-v(n-ninc,3)/v(n-ninc,1)
                r4(m)=v(n,4)/v(n,1)-v(n-ninc,4)/v(n-ninc,1)
                r5(m)=ps(n)-ps(n-ninc)
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1p1,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                if(abs(r1(m-ninc)*r1(m)*r1(m+ninc)).le.tiny(1.)) then
                   r1(m)=1.
                else
                   r1(m)=r1(m)/r1(m-ninc)
                endif
                if(abs(r2(m-ninc)*r2(m)*r2(m+ninc)).le.tiny(1.)) then
                   r2(m)=1.
                else
                   r2(m)=r2(m)/r2(m-ninc)
                endif
                if(abs(r3(m-ninc)*r3(m)*r3(m+ninc)).le.tiny(1.)) then
                   r3(m)=1.
                else
                   r3(m)=r3(m)/r3(m-ninc)
                endif
                if(abs(r4(m-ninc)*r4(m)*r4(m+ninc)).le.tiny(1.)) then
                   r4(m)=1.
                else
                   r4(m)=r4(m)/r4(m-ninc)
                endif
                if(abs(r5(m-ninc)*r5(m)*r5(m+ninc)).le.tiny(1.)) then
                   r5(m)=1.
                else
                   r5(m)=r5(m)/r5(m-ninc)
                endif
             enddo
          enddo
       enddo

       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1p1,j,k)
             ind2 = indc(i2-2,j,k)
             do n=ind1,ind2
                m=n-n0c
                cal=0.5*cvi(n-ninc)*(1.-xk)/(cvi(n-ninc)+cvi(n))
                cbl=0.5*cvi(n-ninc)*(1.+xk*cvi(n-ninc)/cvi(n))/(cvi(n-ninc)+cvi(n))

                rhol(m)=v(n-ninc,1)+(cal*phi(r1(m)   )*(v(n-ninc,1)-v(n-2*ninc,1)) &
                                    +cbl*phi(1./r1(m))*(v(n     ,1)-v(n-ninc  ,1)))
                ul(m)=v(n-ninc,2)/v(n-ninc,1) + ( &
                     cal*phi(r2(m)   )*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
                    +cbl*phi(1./r2(m))*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
                vl(m)=v(n-ninc,3)/v(n-ninc,1) + ( &
                     cal*phi(r3(m)   )*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
                    +cbl*phi(1./r3(m))*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
                wl(m)=v(n-ninc,4)/v(n-ninc,1) + ( &
                     cal*phi(r4(m)   )*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
                    +cbl*phi(1./r4(m))*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
                pl(m)=ps(n-ninc) + (cal*phi(r5(m)   )*(ps(n-ninc)-ps(n-2*ninc)) &
                    +cbl*phi(1./r5(m))*(ps(n     )-ps(n-  ninc)))
!
                car=0.5*cvi(n)*(1.+xk)/(cvi(n)+cvi(n+ninc))
                cbr=0.5*cvi(n)*(1.-xk*cvi(n)/cvi(n+ninc))/(cvi(n)+cvi(n+ninc))

                rhor(m)=v(n,1)-(car*phi(r1(m+ninc   ))*(v(n,1)     -v(n-ninc,1)) &
                               +cbr*phi(1./r1(m+ninc))*(v(n+ninc,1)-v(n     ,1)))
                ur(m)=v(n,2)/v(n,1) - ( &
                     car*phi(r2(m+ninc   ))*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
                    +cbr*phi(1./r2(m+ninc))*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
                vr(m)=v(n,3)/v(n,1) - ( &
                     car*phi(r3(m+ninc   ))*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
                    +cbr*phi(1./r3(m+ninc))*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
                wr(m)=v(n,4)/v(n,1) - ( &
                     car*phi(r4(m+ninc   ))*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
                    +cbr*phi(1./r4(m+ninc))*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
                prr(m)=ps(n) - (car*phi(r5(m+ninc   ))*(ps(n)     -ps(n-ninc)) &
                               +cbr*phi(1./r5(m+ninc))*(ps(n+ninc)-ps(n     )))
             enddo
          enddo
       enddo

       do k=k1,k2m1
          ind1 = indc(i2m1,j1  ,k)
          ind2 = indc(i2m1,j2m1,k)
          do n=ind1,ind2,ncj
             m=n-n0c
             cal=0.5*cvi(n-ninc)*(1.-xk)/(cvi(n-ninc)+cvi(n))
             cbl=0.5*cvi(n-ninc)*(1.+xk*cvi(n-ninc)/cvi(n))/(cvi(n-ninc)+cvi(n))

             rhol(m)=v(n-ninc,1)+muscl*(cal*(v(n-ninc,1)-v(n-2*ninc,1)) &
                                       +cbl*(v(n     ,1)-v(n-ninc  ,1)))
             ul(m)=v(n-ninc,2)/v(n-ninc,1) + muscl*( &
                  cal*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
                 +cbl*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
             vl(m)=v(n-ninc,3)/v(n-ninc,1) + muscl*( &
                  cal*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
                 +cbl*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
             wl(m)=v(n-ninc,4)/v(n-ninc,1) + muscl*( &
                  cal*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
                 +cbl*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
             pl(m)=ps(n-ninc) + muscl*(cal*(ps(n-ninc)-ps(n-2*ninc)) &
                                      +cbl*(ps(n     )-ps(n-  ninc)))
!
             car=0.5*cvi(n)*(1.+xk)/(cvi(n)+cvi(n+ninc))
             cbr=0.5*cvi(n)*(1.-xk*cvi(n)/cvi(n+ninc))/(cvi(n)+cvi(n+ninc))
             rhor(m)=v(n,1)-muscl*(car*(v(n,1)     -v(n-ninc,1)) &
                                  +cbr*(v(n+ninc,1)-v(n     ,1)))
             ur(m)=v(n,2)/v(n,1) - muscl*( &
                  car*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
                 +cbr*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
             vr(m)=v(n,3)/v(n,1) - muscl*( &
                  car*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
                 +cbr*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
             wr(m)=v(n,4)/v(n,1) - muscl*( &
                  car*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
                 +cbr*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
             prr(m)=ps(n) - muscl*(car*(ps(n)     -ps(n-ninc)) &
                                  +cbr*(ps(n+ninc)-ps(n     )))
          enddo
       enddo

    else

       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1p1,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                cal=0.5*cvi(n-ninc)*(1.-xk)/(cvi(n-ninc)+cvi(n))
                cbl=0.5*cvi(n-ninc)*(1.+xk*cvi(n-ninc)/cvi(n))/ (cvi(n-ninc)+cvi(n))

                rhol(m)=v(n-ninc,1)+muscl*(cal*(v(n-ninc,1)-v(n-2*ninc,1)) &
                                          +cbl*(v(n     ,1)-v(n-ninc  ,1)))
                ul(m)=v(n-ninc,2)/v(n-ninc,1) + muscl*( &
                     cal*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
                    +cbl*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
                vl(m)=v(n-ninc,3)/v(n-ninc,1) + muscl*( &
                     cal*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
                    +cbl*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
                wl(m)=v(n-ninc,4)/v(n-ninc,1) + muscl*( &
                     cal*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
                    +cbl*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
                pl(m)=ps(n-ninc) + muscl*(cal*(ps(n-ninc)-ps(n-2*ninc)) &
                                         +cbl*(ps(n     )-ps(n-  ninc)))
!
                car=0.5*cvi(n)*(1.+xk)/(cvi(n)+cvi(n+ninc))
                cbr=0.5*cvi(n)*(1.-xk*cvi(n)/cvi(n+ninc))/(cvi(n)+cvi(n+ninc))

                rhor(m)=v(n,1)-muscl*(car*(v(n,1)     -v(n-ninc,1)) &
                                     +cbr*(v(n+ninc,1)-v(n     ,1)))
                ur(m)=v(n,2)/v(n,1) - muscl*( &
                     car*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
                    +cbr*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
                vr(m)=v(n,3)/v(n,1) - muscl*( &
                     car*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
                    +cbr*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
                wr(m)=v(n,4)/v(n,1) - muscl*( &
                     car*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
                    +cbr*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
                prr(m)=ps(n) - muscl*(car*(ps(n)     -ps(n-ninc)) &
                                     +cbr*(ps(n+ninc)-ps(n     )))
             enddo
          enddo
       enddo
    endif
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1p1,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             cnds=sqrt(sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3))
!        calcul des etats gauche et droit
             al=sqrt(gam*pl(m)/rhol(m))
             ar=sqrt(gam*prr(m)/rhor(m))
             hl=al*al/gam1+0.5*(ul(m)**2+vl(m)**2+wl(m)**2)
             hr=ar*ar/gam1+0.5*(ur(m)**2+vr(m)**2+wr(m)**2)
!        calcul de la vitesse du son a l'interface
             ai=sqrt(al*ar)
!        calcul des nombres de Mach gauche et droit
             cml=(ul(m)*sn(m,kdir,1)+vl(m)*sn(m,kdir,2)+wl(m)*sn(m,kdir,3))/(ai*cnds)
             cmr=(ur(m)*sn(m,kdir,1)+vr(m)*sn(m,kdir,2)+wr(m)*sn(m,kdir,3))/(ai*cnds)
!        calcul de la pression statique a l'interface (x2)
             psi=2.*(pl(m)*fpp(cml)+prr(m)*fpm(cmr))-2.*pinfl
!        calcul du nombre de Mach a l'interface
             cmi=fmp(cml)+fmm(cmr)
!        calcul du flux de masse
             dm=0.5*ai*(cmi*(rhol(m)+rhor(m))-abs(cmi)*(rhor(m)-rhol(m)))*cnds
!        calcul du flux numerique
             hi1=dm
             hi2=dm*(ul(m)+ur(m))-abs(dm)*(ur(m)-ul(m))+psi*sn(m,kdir,1) &
                  - (cmui1(n)*toxx(n)+cmui2(n)*toxx(n-ninc))*sn(m,kdir,1)  &
                  - (cmui1(n)*toxy(n)+cmui2(n)*toxy(n-ninc))*sn(m,kdir,2) &
                  - (cmui1(n)*toxz(n)+cmui2(n)*toxz(n-ninc))*sn(m,kdir,3) 
             hi3=dm*(vl(m)+vr(m))-abs(dm)*(vr(m)-vl(m))+psi*sn(m,kdir,2) &
                  - (cmui1(n)*toxy(n)+cmui2(n)*toxy(n-ninc))*sn(m,kdir,1) &
                  - (cmui1(n)*toyy(n)+cmui2(n)*toyy(n-ninc))*sn(m,kdir,2) &
                  - (cmui1(n)*toyz(n)+cmui2(n)*toyz(n-ninc))*sn(m,kdir,3)
             hi4=dm*(wl(m)+wr(m))-abs(dm)*(wr(m)-wl(m))+psi*sn(m,kdir,3) &
                  - (cmui1(n)*toxz(n)+cmui2(n)*toxz(n-ninc))*sn(m,kdir,1) &
                  - (cmui1(n)*toyz(n)+cmui2(n)*toyz(n-ninc))*sn(m,kdir,2) &
                  - (cmui1(n)*tozz(n)+cmui2(n)*tozz(n-ninc))*sn(m,kdir,3)
             hi5=dm*(hl+hr)-abs(dm)*(hr-hl)  &
                  - (cmui1(n)*d(n,1)+cmui2(n)*d(n-ninc,1))*sn(m,kdir,1) &
                  - (cmui1(n)*d(n,2)+cmui2(n)*d(n-ninc,2))*sn(m,kdir,2) &
                  - (cmui1(n)*d(n,3)+cmui2(n)*d(n-ninc,3))*sn(m,kdir,3)
             u(n,1)=u(n,1)-hi1
             u(n,2)=u(n,2)-0.5*hi2
             u(n,3)=u(n,3)-0.5*hi3
             u(n,4)=u(n,4)-0.5*hi4
             u(n,5)=u(n,5)-0.5*hi5
             u(n-ninc,1)=u(n-ninc,1)+hi1
             u(n-ninc,2)=u(n-ninc,2)+0.5*hi2
             u(n-ninc,3)=u(n-ninc,3)+0.5*hi3
             u(n-ninc,4)=u(n-ninc,4)+0.5*hi4
             u(n-ninc,5)=u(n-ninc,5)+0.5*hi5
          enddo
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i1,j1  ,k)
       ind2 = indc(i1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          n1=n-ninc
          fxx=v(n1,2)*(v(n1,2)/v(n1,1))+ps(n-ninc)-pinfl-toxx(n1)
          fxy=v(n1,3)*(v(n1,2)/v(n1,1))  -toxy(n1)
          fxz=v(n1,4)*(v(n1,2)/v(n1,1))  -toxz(n1)
          fyy=v(n1,3)*(v(n1,3)/v(n1,1))+ps(n-ninc)-pinfl-toyy(n1)
          fyz=v(n1,4)*(v(n1,3)/v(n1,1))  -toyz(n1)
          fzz=v(n1,4)*(v(n1,4)/v(n1,1))+ps(n-ninc)-pinfl-tozz(n1)
          fex=((v(n1,5)+ps(n-ninc)-pinfl-toxx(n1))*v(n1,2) &
               -toxy(n1)*v(n1,3)-toxz(n1)*v(n1,4))/v(n1,1)-qcx(n1)
          fey=((v(n1,5)+ps(n-ninc)-pinfl-toyy(n1))*v(n1,3) &
               -toxy(n1)*v(n1,2)-toyz(n1)*v(n1,4))/v(n1,1)-qcy(n1)
          fez=((v(n1,5)+ps(n-ninc)-pinfl-tozz(n1))*v(n1,4) &
               -toxz(n1)*v(n1,2)-toyz(n1)*v(n1,3))/v(n1,1)-qcz(n1)
!
          si0= v(n-ninc,2)*sn(m,kdir,1) &
              +v(n-ninc,3)*sn(m,kdir,2) &
              +v(n-ninc,4)*sn(m,kdir,3)
          si1= fxx*sn(m,kdir,1) &
              +fxy*sn(m,kdir,2) &
              +fxz*sn(m,kdir,3)
          si2= fxy*sn(m,kdir,1) &
              +fyy*sn(m,kdir,2) &
              +fyz*sn(m,kdir,3)
          si3= fxz*sn(m,kdir,1) &
              +fyz*sn(m,kdir,2) &
              +fzz*sn(m,kdir,3)
          si4= fex*sn(m,kdir,1) &
              +fey*sn(m,kdir,2) &
              +fez*sn(m,kdir,3)
          u(n,1)=u(n,1)-si0
          u(n,2)=u(n,2)-si1
          u(n,3)=u(n,3)-si2
          u(n,4)=u(n,4)-si3
          u(n,5)=u(n,5)-si4
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i2,j1  ,k)
       ind2 = indc(i2,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          fxx=v(n,2)*(v(n,2)/v(n,1))+ps(n)-pinfl-toxx(n)
          fxy=v(n,3)*(v(n,2)/v(n,1))  -toxy(n)
          fxz=v(n,4)*(v(n,2)/v(n,1))  -toxz(n)
          fyy=v(n,3)*(v(n,3)/v(n,1))+ps(n)-pinfl-toyy(n)
          fyz=v(n,4)*(v(n,3)/v(n,1))  -toyz(n)
          fzz=v(n,4)*(v(n,4)/v(n,1))+ps(n)-pinfl-tozz(n)
          fex=((v(n,5)+ps(n)-pinfl-toxx(n))*v(n,2) &
               -toxy(n)*v(n,3)-toxz(n)*v(n,4))/v(n,1)-qcx(n)
          fey=((v(n,5)+ps(n)-pinfl-toyy(n))*v(n,3) &
               -toxy(n)*v(n,2)-toyz(n)*v(n,4))/v(n,1)-qcy(n)
          fez=((v(n,5)+ps(n)-pinfl-tozz(n))*v(n,4) &
               -toxz(n)*v(n,2)-toyz(n)*v(n,3))/v(n,1)-qcz(n)
!
          si0= v(n,2)*sn(m,kdir,1) &
               +v(n,3)*sn(m,kdir,2) &
               +v(n,4)*sn(m,kdir,3)
          si1= fxx*sn(m,kdir,1) &
               +fxy*sn(m,kdir,2) &
               +fxz*sn(m,kdir,3)
          si2= fxy*sn(m,kdir,1) &
               +fyy*sn(m,kdir,2) &
               +fyz*sn(m,kdir,3)
          si3= fxz*sn(m,kdir,1) &
               +fyz*sn(m,kdir,2) &
               +fzz*sn(m,kdir,3)
          si4= fex*sn(m,kdir,1) &
               +fey*sn(m,kdir,2) &
               +fez*sn(m,kdir,3)
          u(n-ninc,1)=u(n-ninc,1)+si0
          u(n-ninc,2)=u(n-ninc,2)+si1
          u(n-ninc,3)=u(n-ninc,3)+si2
          u(n-ninc,4)=u(n-ninc,4)+si3
          u(n-ninc,5)=u(n-ninc,5)+si4
       enddo
    enddo
!
!------direction j----------------------------------------------
!
    kdir=2
    ninc=ncj
!
!-----definition des variables extrapolees------------------------
!
    if(ilim.eq.1) then
!
       do k=k1,k2m1
          do j=j1,j2
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                r1(m)=v(n,1)-v(n-ninc,1)
                r2(m)=v(n,2)/v(n,1)-v(n-ninc,2)/v(n-ninc,1)
                r3(m)=v(n,3)/v(n,1)-v(n-ninc,3)/v(n-ninc,1)
                r4(m)=v(n,4)/v(n,1)-v(n-ninc,4)/v(n-ninc,1)
                r5(m)=ps(n)-ps(n-ninc)
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          do j=j1p1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                if(abs(r1(m-ninc)*r1(m)*r1(m+ninc)).le.tiny(1.)) then
                   r1(m)=1.
                else
                   r1(m)=r1(m)/r1(m-ninc)
                endif
                if(abs(r2(m-ninc)*r2(m)*r2(m+ninc)).le.tiny(1.)) then
                   r2(m)=1.
                else
                   r2(m)=r2(m)/r2(m-ninc)
                endif
                if(abs(r3(m-ninc)*r3(m)*r3(m+ninc)).le.tiny(1.)) then
                   r3(m)=1.
                else
                   r3(m)=r3(m)/r3(m-ninc)
                endif
                if(abs(r4(m-ninc)*r4(m)*r4(m+ninc)).le.tiny(1.)) then
                   r4(m)=1.
                else
                   r4(m)=r4(m)/r4(m-ninc)
                endif
                if(abs(r5(m-ninc)*r5(m)*r5(m+ninc)).le.tiny(1.)) then
                   r5(m)=1.
                else
                   r5(m)=r5(m)/r5(m-ninc)
                endif
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          do j=j1p1,j2-2
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                cal=0.5*cvj(n-ninc)*(1.-xk)/(cvj(n-ninc)+cvj(n))
                cbl=0.5*cvj(n-ninc)*(1.+xk*cvj(n-ninc)/cvj(n))/(cvj(n-ninc)+cvj(n))

                rhol(m)=v(n-ninc,1)+(cal*phi(r1(m)   )*(v(n-ninc,1)-v(n-2*ninc,1)) &
                     +cbl*phi(1./r1(m))*(v(n     ,1)-v(n-ninc  ,1)))
                ul(m)=v(n-ninc,2)/v(n-ninc,1) + ( &
                     cal*phi(r2(m)   )*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
                     +cbl*phi(1./r2(m))*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
                vl(m)=v(n-ninc,3)/v(n-ninc,1) + ( &
                     cal*phi(r3(m)   )*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
                     +cbl*phi(1./r3(m))*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
                wl(m)=v(n-ninc,4)/v(n-ninc,1) + ( &
                     cal*phi(r4(m)   )*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
                     +cbl*phi(1./r4(m))*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
                pl(m)=ps(n-ninc) + (cal*phi(r5(m)   )*(ps(n-ninc)-ps(n-2*ninc)) &
                     +cbl*phi(1./r5(m))*(ps(n     )-ps(n-  ninc)))
!
                car=0.5*cvj(n)*(1.+xk)/(cvj(n)+cvj(n+ninc))
                cbr=0.5*cvj(n)*(1.-xk*cvj(n)/cvj(n+ninc))/(cvj(n)+cvj(n+ninc))

                rhor(m)=v(n,1)-(car*phi(r1(m+ninc   ))*(v(n,1)     -v(n-ninc,1)) &
                     +cbr*phi(1./r1(m+ninc))*(v(n+ninc,1)-v(n     ,1)))
                ur(m)=v(n,2)/v(n,1) - ( &
                     car*phi(r2(m+ninc   ))*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
                     +cbr*phi(1./r2(m+ninc))*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
                vr(m)=v(n,3)/v(n,1) - ( &
                     car*phi(r3(m+ninc   ))*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
                     +cbr*phi(1./r3(m+ninc))*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
                wr(m)=v(n,4)/v(n,1) - ( &
                     car*phi(r4(m+ninc   ))*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
                     +cbr*phi(1./r4(m+ninc))*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
                prr(m)=ps(n) - (car*phi(r5(m+ninc   ))*(ps(n)     -ps(n-ninc)) &
                     +cbr*phi(1./r5(m+ninc))*(ps(n+ninc)-ps(n     )))
             enddo
          enddo
       enddo

       do k=k1,k2m1
          ind1 = indc(i1  ,j2m1,k)
          ind2 = indc(i2m1,j2m1,k)
          do n=ind1,ind2
             m=n-n0c
             cal=0.5*cvj(n-ninc)*(1.-xk)/(cvj(n-ninc)+cvj(n))
             cbl=0.5*cvj(n-ninc)*(1.+xk*cvj(n-ninc)/cvj(n))/(cvj(n-ninc)+cvj(n))

             rhol(m)=v(n-ninc,1)+muscl*(cal*(v(n-ninc,1)-v(n-2*ninc,1)) &
                  +cbl*(v(n     ,1)-v(n-ninc  ,1)))
             ul(m)=v(n-ninc,2)/v(n-ninc,1) + muscl*( &
                  cal*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
                  +cbl*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
             vl(m)=v(n-ninc,3)/v(n-ninc,1) + muscl*( &
                  cal*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
                  +cbl*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
             wl(m)=v(n-ninc,4)/v(n-ninc,1) + muscl*( &
                  cal*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
                  +cbl*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
             pl(m)=ps(n-ninc) + muscl*(cal*(ps(n-ninc)-ps(n-2*ninc)) &
                  +cbl*(ps(n     )-ps(n-  ninc)))
!
             car=0.5*cvj(n)*(1.+xk)/(cvj(n)+cvj(n+ninc))
             cbr=0.5*cvj(n)*(1.-xk*cvj(n)/cvj(n+ninc))/(cvj(n)+cvj(n+ninc))

             rhor(m)=v(n,1)-muscl*(car*(v(n,1)     -v(n-ninc,1)) &
                  +cbr*(v(n+ninc,1)-v(n     ,1)))
             ur(m)=v(n,2)/v(n,1) - muscl*( &
                  car*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
                  +cbr*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
             vr(m)=v(n,3)/v(n,1) - muscl*( &
                  car*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
                  +cbr*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
             wr(m)=v(n,4)/v(n,1) - muscl*( &
                  car*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
                  +cbr*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
             prr(m)=ps(n) - muscl*(car*(ps(n)     -ps(n-ninc)) &
                  +cbr*(ps(n+ninc)-ps(n     )))
          enddo
       enddo
!
    else
!
       do k=k1,k2m1
          do j=j1p1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                cal=0.5*cvj(n-ninc)*(1.-xk)/(cvj(n-ninc)+cvj(n))
                cbl=0.5*cvj(n-ninc)*(1.+xk*cvj(n-ninc)/cvj(n))/(cvj(n-ninc)+cvj(n))

                rhol(m)=v(n-ninc,1)+muscl*(cal*(v(n-ninc,1)-v(n-2*ninc,1)) &
                     +cbl*(v(n     ,1)-v(n-ninc  ,1)))
                ul(m)=v(n-ninc,2)/v(n-ninc,1) + muscl*( &
                     cal*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
                     +cbl*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
                vl(m)=v(n-ninc,3)/v(n-ninc,1) + muscl*( &
                     cal*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
                     +cbl*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
                wl(m)=v(n-ninc,4)/v(n-ninc,1) + muscl*( &
                     cal*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
                     +cbl*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
                pl(m)=ps(n-ninc) + muscl*(cal*(ps(n-ninc)-ps(n-2*ninc)) &
                     +cbl*(ps(n     )-ps(n-  ninc)))
!
                car=0.5*cvj(n)*(1.+xk)/(cvj(n)+cvj(n+ninc))
                cbr=0.5*cvj(n)*(1.-xk*cvj(n)/cvj(n+ninc))/(cvj(n)+cvj(n+ninc))

                rhor(m)=v(n,1)-muscl*(car*(v(n,1)     -v(n-ninc,1)) &
                     +cbr*(v(n+ninc,1)-v(n     ,1)))
                ur(m)=v(n,2)/v(n,1) - muscl*( &
                     car*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
                     +cbr*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
                vr(m)=v(n,3)/v(n,1) - muscl*( &
                     car*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
                     +cbr*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
                wr(m)=v(n,4)/v(n,1) - muscl*( &
                     car*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
                     +cbr*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
                prr(m)=ps(n) - muscl*(car*(ps(n)     -ps(n-ninc)) &
                     +cbr*(ps(n+ninc)-ps(n     )))
             enddo
          enddo
       enddo
    endif
!
    do k=k1,k2m1
       do j=j1p1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             cnds=sqrt(sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3))
!        calcul des etats gauche et droit
             al=sqrt(gam*pl(m)/rhol(m))
             ar=sqrt(gam*prr(m)/rhor(m))
             hl=al*al/gam1+0.5*(ul(m)**2+vl(m)**2+wl(m)**2)
             hr=ar*ar/gam1+0.5*(ur(m)**2+vr(m)**2+wr(m)**2)
!        calcul de la vitesse du son a l'interface
             ai=sqrt(al*ar)
!        calcul des nombres de Mach gauche et droit
             cml=(ul(m)*sn(m,kdir,1)+vl(m)*sn(m,kdir,2)+wl(m)*sn(m,kdir,3))/(al*cnds)
             cmr=(ur(m)*sn(m,kdir,1)+vr(m)*sn(m,kdir,2)+wr(m)*sn(m,kdir,3))/(ar*cnds)
!        calcul de la pression statique a l'interface (x2)
             psi=2.*(pl(m)*fpp(cml)+prr(m)*fpm(cmr))-2.*pinfl
!        calcul du nombre de Mach a l'interface
             cmi=fmp(cml)+fmm(cmr)
!        calcul du flux de masse
             dm=0.5*ai*(cmi*(rhol(m)+rhor(m))-abs(cmi)*(rhor(m)-rhol(m)))*cnds
!        calcul du flux numerique
             hj1=dm
             hj2=dm*(ul(m)+ur(m))-abs(dm)*(ur(m)-ul(m))+psi*sn(m,kdir,1) &
                  - (cmuj1(n)*toxx(n)+cmuj2(n)*toxx(n-ninc))*sn(m,kdir,1)  &
                  - (cmuj1(n)*toxy(n)+cmuj2(n)*toxy(n-ninc))*sn(m,kdir,2) &
                  - (cmuj1(n)*toxz(n)+cmuj2(n)*toxz(n-ninc))*sn(m,kdir,3)
             hj3=dm*(vl(m)+vr(m))-abs(dm)*(vr(m)-vl(m))+psi*sn(m,kdir,2) &
                  - (cmuj1(n)*toxy(n)+cmuj2(n)*toxy(n-ninc))*sn(m,kdir,1) &
                  - (cmuj1(n)*toyy(n)+cmuj2(n)*toyy(n-ninc))*sn(m,kdir,2) &
                  - (cmuj1(n)*toyz(n)+cmuj2(n)*toyz(n-ninc))*sn(m,kdir,3)
             hj4=dm*(wl(m)+wr(m))-abs(dm)*(wr(m)-wl(m))+psi*sn(m,kdir,3) &
                  - (cmuj1(n)*toxz(n)+cmuj2(n)*toxz(n-ninc))*sn(m,kdir,1) &
                  - (cmuj1(n)*toyz(n)+cmuj2(n)*toyz(n-ninc))*sn(m,kdir,2) &
                  - (cmuj1(n)*tozz(n)+cmuj2(n)*tozz(n-ninc))*sn(m,kdir,3)
             hj5=dm*(hl+hr)-abs(dm)*(hr-hl)  &
                  - (cmuj1(n)*d(n,1)+cmuj2(n)*d(n-ninc,1))*sn(m,kdir,1) &
                  - (cmuj1(n)*d(n,2)+cmuj2(n)*d(n-ninc,2))*sn(m,kdir,2) &
                  - (cmuj1(n)*d(n,3)+cmuj2(n)*d(n-ninc,3))*sn(m,kdir,3)
             u(n,1)=u(n,1)-hj1
             u(n,2)=u(n,2)-0.5*hj2
             u(n,3)=u(n,3)-0.5*hj3
             u(n,4)=u(n,4)-0.5*hj4
             u(n,5)=u(n,5)-0.5*hj5
             u(n-ninc,1)=u(n-ninc,1)+hj1
             u(n-ninc,2)=u(n-ninc,2)+0.5*hj2
             u(n-ninc,3)=u(n-ninc,3)+0.5*hj3
             u(n-ninc,4)=u(n-ninc,4)+0.5*hj4
             u(n-ninc,5)=u(n-ninc,5)+0.5*hj5
          enddo
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j1,k)
       ind2 = indc(i2m1,j1,k)
       do n=ind1,ind2
          m=n-n0c
          n1=n-ninc
          fxx=v(n1,2)*(v(n1,2)/v(n1,1))+ps(n-ninc)-pinfl-toxx(n1)
          fxy=v(n1,3)*(v(n1,2)/v(n1,1))  -toxy(n1)
          fxz=v(n1,4)*(v(n1,2)/v(n1,1))  -toxz(n1)
          fyy=v(n1,3)*(v(n1,3)/v(n1,1))+ps(n-ninc)-pinfl-toyy(n1)
          fyz=v(n1,4)*(v(n1,3)/v(n1,1))  -toyz(n1)
          fzz=v(n1,4)*(v(n1,4)/v(n1,1))+ps(n-ninc)-pinfl-tozz(n1)
          fex=((v(n1,5)+ps(n-ninc)-pinfl-toxx(n1))*v(n1,2) &
               -toxy(n1)*v(n1,3)-toxz(n1)*v(n1,4))/v(n1,1)-qcx(n1)
          fey=((v(n1,5)+ps(n-ninc)-pinfl-toyy(n1))*v(n1,3) &
               -toxy(n1)*v(n1,2)-toyz(n1)*v(n1,4))/v(n1,1)-qcy(n1)
          fez=((v(n1,5)+ps(n-ninc)-pinfl-tozz(n1))*v(n1,4) &
               -toxz(n1)*v(n1,2)-toyz(n1)*v(n1,3))/v(n1,1)-qcz(n1)
!
          sj0= v(n-ninc,2)*sn(m,kdir,1) &
              +v(n-ninc,3)*sn(m,kdir,2) &
              +v(n-ninc,4)*sn(m,kdir,3)
          sj1= fxx*sn(m,kdir,1) &
              +fxy*sn(m,kdir,2) &
              +fxz*sn(m,kdir,3)
          sj2= fxy*sn(m,kdir,1) &
              +fyy*sn(m,kdir,2) &
              +fyz*sn(m,kdir,3)
          sj3= fxz*sn(m,kdir,1) &
              +fyz*sn(m,kdir,2) &
              +fzz*sn(m,kdir,3)
          sj4= fex*sn(m,kdir,1) &
              +fey*sn(m,kdir,2) &
              +fez*sn(m,kdir,3)
          u(n,1)=u(n,1)-sj0
          u(n,2)=u(n,2)-sj1
          u(n,3)=u(n,3)-sj2
          u(n,4)=u(n,4)-sj3
          u(n,5)=u(n,5)-sj4
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j2,k)
       ind2 = indc(i2m1,j2,k)
       do n=ind1,ind2
          m=n-n0c
          fxx=v(n,2)*(v(n,2)/v(n,1))+ps(n)-pinfl-toxx(n)
          fxy=v(n,3)*(v(n,2)/v(n,1))  -toxy(n)
          fxz=v(n,4)*(v(n,2)/v(n,1))  -toxz(n)
          fyy=v(n,3)*(v(n,3)/v(n,1))+ps(n)-pinfl-toyy(n)
          fyz=v(n,4)*(v(n,3)/v(n,1))  -toyz(n)
          fzz=v(n,4)*(v(n,4)/v(n,1))+ps(n)-pinfl-tozz(n)
          fex=((v(n,5)+ps(n)-pinfl-toxx(n))*v(n,2) &
               -toxy(n)*v(n,3)-toxz(n)*v(n,4))/v(n,1)-qcx(n)
          fey=((v(n,5)+ps(n)-pinfl-toyy(n))*v(n,3) &
               -toxy(n)*v(n,2)-toyz(n)*v(n,4))/v(n,1)-qcy(n)
          fez=((v(n,5)+ps(n)-pinfl-tozz(n))*v(n,4) &
               -toxz(n)*v(n,2)-toyz(n)*v(n,3))/v(n,1)-qcz(n)
!
          sj0= v(n,2)*sn(m,kdir,1) &
              +v(n,3)*sn(m,kdir,2) &
              +v(n,4)*sn(m,kdir,3)
          sj1= fxx*sn(m,kdir,1) &
              +fxy*sn(m,kdir,2) &
              +fxz*sn(m,kdir,3)
          sj2= fxy*sn(m,kdir,1) &
              +fyy*sn(m,kdir,2) &
              +fyz*sn(m,kdir,3)
          sj3= fxz*sn(m,kdir,1) &
              +fyz*sn(m,kdir,2) &
              +fzz*sn(m,kdir,3)
          sj4= fex*sn(m,kdir,1) &
              +fey*sn(m,kdir,2) &
              +fez*sn(m,kdir,3)
          u(n-ninc,1)=u(n-ninc,1)+sj0
          u(n-ninc,2)=u(n-ninc,2)+sj1
          u(n-ninc,3)=u(n-ninc,3)+sj2
          u(n-ninc,4)=u(n-ninc,4)+sj3
          u(n-ninc,5)=u(n-ninc,5)+sj4
       enddo
    enddo
!
!------direction k----------------------------------------------
!
    if(equat(3:4).eq.'3d') then
       kdir=3
       ninc=nck
!
       if(ilim.eq.1) then
!
          do k=k1,k2
             do j=j1,j2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   r1(m)=v(n,1)-v(n-ninc,1)
                   r2(m)=v(n,2)/v(n,1)-v(n-ninc,2)/v(n-ninc,1)
                   r3(m)=v(n,3)/v(n,1)-v(n-ninc,3)/v(n-ninc,1)
                   r4(m)=v(n,4)/v(n,1)-v(n-ninc,4)/v(n-ninc,1)
                   r5(m)=ps(n)-ps(n-ninc)
                enddo
             enddo
          enddo
!
          do k=k1p1,k2m1
             do j=j1,j2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   if(abs(r1(m-ninc)*r1(m)*r1(m+ninc)).le.tiny(1.)) then
                      r1(m)=1.
                   else
                      r1(m)=r1(m)/r1(m-ninc)
                   endif
                   if(abs(r2(m-ninc)*r2(m)*r2(m+ninc)).le.tiny(1.)) then
                      r2(m)=1.
                   else
                      r2(m)=r2(m)/r2(m-ninc)
                   endif
                   if(abs(r3(m-ninc)*r3(m)*r3(m+ninc)).le.tiny(1.)) then
                      r3(m)=1.
                   else
                      r3(m)=r3(m)/r3(m-ninc)
                   endif
                   if(abs(r4(m-ninc)*r4(m)*r4(m+ninc)).le.tiny(1.)) then
                      r4(m)=1.
                   else
                      r4(m)=r4(m)/r4(m-ninc)
                   endif
                   if(abs(r5(m-ninc)*r5(m)*r5(m+ninc)).le.tiny(1.)) then
                      r5(m)=1.
                   else
                      r5(m)=r5(m)/r5(m-ninc)
                   endif
                enddo
             enddo
          enddo

          do k=k1p1,k2-2
             do j=j1,j2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   cal=0.5*cvk(n-ninc)*(1.-xk)/(cvk(n-ninc)+cvk(n))
                   cbl=0.5*cvk(n-ninc)*(1.+xk*cvk(n-ninc)/cvk(n))/(cvk(n-ninc)+cvk(n))

                   rhol(m)=v(n-ninc,1)+(cal*phi(r1(m)   )*(v(n-ninc,1)-v(n-2*ninc,1)) &
                        +cbl*phi(1./r1(m))*(v(n     ,1)-v(n-ninc  ,1)))
                   ul(m)=v(n-ninc,2)/v(n-ninc,1) + ( &
                        cal*phi(r2(m)   )*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
                        +cbl*phi(1./r2(m))*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
                   vl(m)=v(n-ninc,3)/v(n-ninc,1) + ( &
                        cal*phi(r3(m)   )*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
                        +cbl*phi(1./r3(m))*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
                   wl(m)=v(n-ninc,4)/v(n-ninc,1) + ( &
                        cal*phi(r4(m)   )*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
                        +cbl*phi(1./r4(m))*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
                   pl(m)=ps(n-ninc) + (cal*phi(r5(m)   )*(ps(n-ninc)-ps(n-2*ninc)) &
                        +cbl*phi(1./r5(m))*(ps(n     )-ps(n-  ninc)))
!
                   car=0.5*cvk(n)*(1.+xk)/(cvk(n)+cvk(n+ninc))
                   cbr=0.5*cvk(n)*(1.-xk*cvk(n)/cvk(n+ninc))/(cvk(n)+cvk(n+ninc))

                   rhor(m)=v(n,1)-(car*phi(r1(m+ninc   ))*(v(n,1)     -v(n-ninc,1)) &
                        +cbr*phi(1./r1(m+ninc))*(v(n+ninc,1)-v(n     ,1)))
                   ur(m)=v(n,2)/v(n,1) - ( &
                        car*phi(r2(m+ninc   ))*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
                        +cbr*phi(1./r2(m+ninc))*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
                   vr(m)=v(n,3)/v(n,1) - ( &
                        car*phi(r3(m+ninc   ))*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
                        +cbr*phi(1./r3(m+ninc))*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
                   wr(m)=v(n,4)/v(n,1) - ( &
                        car*phi(r4(m+ninc   ))*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
                        +cbr*phi(1./r4(m+ninc))*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
                   prr(m)=ps(n) - (car*phi(r5(m+ninc   ))*(ps(n)     -ps(n-ninc)) &
                        +cbr*phi(1./r5(m+ninc))*(ps(n+ninc)-ps(n     )))
                enddo
             enddo
          enddo

          do k=k1,k2m1
             ind1 = indc(i1  ,j1  ,k2m1)
             ind2 = indc(i2m1,j2m1,k2m1)
             do n=ind1,ind2
                m=n-n0c
                cal=0.5*cvk(n-ninc)*(1.-xk)/(cvk(n-ninc)+cvk(n))
                cbl=0.5*cvk(n-ninc)*(1.+xk*cvk(n-ninc)/cvk(n))/(cvk(n-ninc)+cvk(n))

                rhol(m)=v(n-ninc,1)+muscl*(cal*(v(n-ninc,1)-v(n-2*ninc,1)) &
                     +cbl*(v(n     ,1)-v(n-ninc  ,1)))
                ul(m)=v(n-ninc,2)/v(n-ninc,1) + muscl*( &
                     cal*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
                     +cbl*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
                vl(m)=v(n-ninc,3)/v(n-ninc,1) + muscl*( &
                     cal*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
                     +cbl*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
                wl(m)=v(n-ninc,4)/v(n-ninc,1) + muscl*( &
                     cal*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
                     +cbl*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
                pl(m)=ps(n-ninc) + muscl*(cal*(ps(n-ninc)-ps(n-2*ninc)) &
                     +cbl*(ps(n     )-ps(n-  ninc)))
!
                car=0.5*cvk(n)*(1.+xk)/(cvk(n)+cvk(n+ninc))
                cbr=0.5*cvk(n)*(1.-xk*cvk(n)/cvk(n+ninc))/(cvk(n)+cvk(n+ninc))

                rhor(m)=v(n,1)-muscl*(car*(v(n,1)     -v(n-ninc,1)) &
                     +cbr*(v(n+ninc,1)-v(n     ,1)))
                ur(m)=v(n,2)/v(n,1) - muscl*( &
                     car*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
                     +cbr*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
                vr(m)=v(n,3)/v(n,1) - muscl*( &
                     car*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
                     +cbr*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
                wr(m)=v(n,4)/v(n,1) - muscl*( &
                     car*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
                     +cbr*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
                prr(m)=ps(n) - muscl*(car*(ps(n)     -ps(n-ninc)) &
                     +cbr*(ps(n+ninc)-ps(n     )))
             enddo
          enddo
!
       else
!
          do k=k1p1,k2m1
             do j=j1,j2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   cal=0.5*cvk(n-ninc)*(1.-xk)/(cvk(n-ninc)+cvk(n))
                   cbl=0.5*cvk(n-ninc)*(1.+xk*cvk(n-ninc)/cvk(n))/(cvk(n-ninc)+cvk(n))

                   rhol(m)=v(n-ninc,1)+muscl*(cal*(v(n-ninc,1)-v(n-2*ninc,1)) &
                        +cbl*(v(n     ,1)-v(n-ninc  ,1)))
                   ul(m)=v(n-ninc,2)/v(n-ninc,1) + muscl*( &
                        cal*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
                        +cbl*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
                   vl(m)=v(n-ninc,3)/v(n-ninc,1) + muscl*( &
                        cal*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
                        +cbl*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
                   wl(m)=v(n-ninc,4)/v(n-ninc,1) + muscl*( &
                        cal*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
                        +cbl*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
                   pl(m)=ps(n-ninc) + muscl*(cal*(ps(n-ninc)-ps(n-2*ninc)) &
                        +cbl*(ps(n     )-ps(n-  ninc)))
!
                   car=0.5*cvk(n)*(1.+xk)/(cvk(n)+cvk(n+ninc))
                   cbr=0.5*cvk(n)*(1.-xk*cvk(n)/cvk(n+ninc))/(cvk(n)+cvk(n+ninc))

                   rhor(m)=v(n,1)-muscl*(car*(v(n,1)     -v(n-ninc,1)) &
                        +cbr*(v(n+ninc,1)-v(n     ,1)))
                   ur(m)=v(n,2)/v(n,1) - muscl*( &
                        car*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
                        +cbr*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
                   vr(m)=v(n,3)/v(n,1) - muscl*( &
                        car*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
                        +cbr*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
                   wr(m)=v(n,4)/v(n,1) - muscl*( &
                        car*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
                        +cbr*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
                   prr(m)=ps(n) - muscl*(car*(ps(n)     -ps(n-ninc)) &
                        +cbr*(ps(n+ninc)-ps(n     )))
                enddo
             enddo
          enddo
       endif
!
       do k=k1p1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                cnds=sqrt(sn(m,kdir,1)*sn(m,kdir,1)+ &
                     sn(m,kdir,2)*sn(m,kdir,2)+ &
                     sn(m,kdir,3)*sn(m,kdir,3))
!        calcul des etats gauche et droit
                al=sqrt(gam*pl(m)/rhol(m))
                ar=sqrt(gam*prr(m)/rhor(m))
                hl=al*al/gam1+0.5*(ul(m)**2+vl(m)**2+wl(m)**2)
                hr=ar*ar/gam1+0.5*(ur(m)**2+vr(m)**2+wr(m)**2)
!        calcul de la vitesse du son a l'interface
                ai=sqrt(al*ar)
!        calcul des nombres de Mach gauche et droit
                cml=(ul(m)*sn(m,kdir,1)+vl(m)*sn(m,kdir,2)+wl(m)*sn(m,kdir,3))/(al*cnds)
                cmr=(ur(m)*sn(m,kdir,1)+vr(m)*sn(m,kdir,2)+wr(m)*sn(m,kdir,3))/(ar*cnds)
!        calcul de la pression statique a l'interface (x2)
                psi=2.*(pl(m)*fpp(cml)+prr(m)*fpm(cmr))-2.*pinfl
!        calcul du nombre de Mach a l'interface
                cmi=fmp(cml)+fmm(cmr)
!        calcul du flux de masse
                dm=0.5*ai*(cmi*(rhol(m)+rhor(m))-abs(cmi)*(rhor(m)-rhol(m)))*cnds
!        calcul du flux numerique
                hk1=dm
                hk2=dm*(ul(m)+ur(m))-abs(dm)*(ur(m)-ul(m))+psi*sn(m,kdir,1) &
                     - (cmuk1(n)*toxx(n)+cmuk2(n)*toxx(n-ninc))*sn(m,kdir,1) &
                     - (cmuk1(n)*toxy(n)+cmuk2(n)*toxy(n-ninc))*sn(m,kdir,2) &
                     - (cmuk1(n)*toxz(n)+cmuk2(n)*toxz(n-ninc))*sn(m,kdir,3)
                hk3=dm*(vl(m)+vr(m))-abs(dm)*(vr(m)-vl(m))+psi*sn(m,kdir,2) &
                     - (cmuk1(n)*toxy(n)+cmuk2(n)*toxy(n-ninc))*sn(m,kdir,1) &
                     - (cmuk1(n)*toyy(n)+cmuk2(n)*toyy(n-ninc))*sn(m,kdir,2) &
                     - (cmuk1(n)*toyz(n)+cmuk2(n)*toyz(n-ninc))*sn(m,kdir,3)
                hk4=dm*(wl(m)+wr(m))-abs(dm)*(wr(m)-wl(m))+psi*sn(m,kdir,3) &
                     - (cmuk1(n)*toxz(n)+cmuk2(n)*toxz(n-ninc))*sn(m,kdir,1) &
                     - (cmuk1(n)*toyz(n)+cmuk2(n)*toyz(n-ninc))*sn(m,kdir,2) &
                     - (cmuk1(n)*tozz(n)+cmuk2(n)*tozz(n-ninc))*sn(m,kdir,3)
                hk5=dm*(hl+hr)-abs(dm)*(hr-hl)  &
                     - (cmuk1(n)*d(n,1)+cmuk2(n)*d(n-ninc,1))*sn(m,kdir,1) &
                     - (cmuk1(n)*d(n,2)+cmuk2(n)*d(n-ninc,2))*sn(m,kdir,2) &
                     - (cmuk1(n)*d(n,3)+cmuk2(n)*d(n-ninc,3))*sn(m,kdir,3)
                u(n,1)=u(n,1)-hk1
                u(n,2)=u(n,2)-0.5*hk2
                u(n,3)=u(n,3)-0.5*hk3
                u(n,4)=u(n,4)-0.5*hk4
                u(n,5)=u(n,5)-0.5*hk5
                u(n-ninc,1)=u(n-ninc,1)+hk1
                u(n-ninc,2)=u(n-ninc,2)+0.5*hk2
                u(n-ninc,3)=u(n-ninc,3)+0.5*hk3
                u(n-ninc,4)=u(n-ninc,4)+0.5*hk4
                u(n-ninc,5)=u(n-ninc,5)+0.5*hk5
             enddo
          enddo
       enddo
!
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k1)
          ind2 = indc(i2m1,j,k1)
          do n=ind1,ind2
             m=n-n0c
             n1=n-ninc
             fxx=v(n1,2)*(v(n1,2)/v(n1,1))+ps(n-ninc)-pinfl-toxx(n1)
             fxy=v(n1,3)*(v(n1,2)/v(n1,1))  -toxy(n1)
             fxz=v(n1,4)*(v(n1,2)/v(n1,1))  -toxz(n1)
             fyy=v(n1,3)*(v(n1,3)/v(n1,1))+ps(n-ninc)-pinfl-toyy(n1)
             fyz=v(n1,4)*(v(n1,3)/v(n1,1))  -toyz(n1)
             fzz=v(n1,4)*(v(n1,4)/v(n1,1))+ps(n-ninc)-pinfl-tozz(n1)
             fex=((v(n1,5)+ps(n-ninc)-pinfl-toxx(n1))*v(n1,2) &
                  -toxy(n1)*v(n1,3)-toxz(n1)*v(n1,4))/v(n1,1)-qcx(n1)
             fey=((v(n1,5)+ps(n-ninc)-pinfl-toyy(n1))*v(n1,3) &
                  -toxy(n1)*v(n1,2)-toyz(n1)*v(n1,4))/v(n1,1)-qcy(n1)
             fez=((v(n1,5)+ps(n-ninc)-pinfl-tozz(n1))*v(n1,4) &
                  -toxz(n1)*v(n1,2)-toyz(n1)*v(n1,3))/v(n1,1)-qcz(n1)
!
             sk0= v(n-ninc,2)*sn(m,kdir,1) &
                 +v(n-ninc,3)*sn(m,kdir,2) &
                 +v(n-ninc,4)*sn(m,kdir,3)
             sk1= fxx*sn(m,kdir,1) &
                 +fxy*sn(m,kdir,2) &
                 +fxz*sn(m,kdir,3)
             sk2= fxy*sn(m,kdir,1) &
                 +fyy*sn(m,kdir,2) &
                 +fyz*sn(m,kdir,3)
             sk3= fxz*sn(m,kdir,1) &
                 +fyz*sn(m,kdir,2) &
                 +fzz*sn(m,kdir,3)
             sk4= fex*sn(m,kdir,1) &
                 +fey*sn(m,kdir,2) &
                 +fez*sn(m,kdir,3)
             u(n,1)=u(n,1)-sk0
             u(n,2)=u(n,2)-sk1
             u(n,3)=u(n,3)-sk2
             u(n,4)=u(n,4)-sk3
             u(n,5)=u(n,5)-sk4
          enddo
       enddo
!
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k2)
          ind2 = indc(i2m1,j,k2)
          do n=ind1,ind2
             m=n-n0c
             fxx=v(n,2)*(v(n,2)/v(n,1))+ps(n)-pinfl-toxx(n)
             fxy=v(n,3)*(v(n,2)/v(n,1))  -toxy(n)
             fxz=v(n,4)*(v(n,2)/v(n,1))  -toxz(n)
             fyy=v(n,3)*(v(n,3)/v(n,1))+ps(n)-pinfl-toyy(n)
             fyz=v(n,4)*(v(n,3)/v(n,1))  -toyz(n)
             fzz=v(n,4)*(v(n,4)/v(n,1))+ps(n)-pinfl-tozz(n)
             fex=((v(n,5)+ps(n)-pinfl-toxx(n))*v(n,2) &
                  -toxy(n)*v(n,3)-toxz(n)*v(n,4))/v(n,1)-qcx(n)
             fey=((v(n,5)+ps(n)-pinfl-toyy(n))*v(n,3) &
                  -toxy(n)*v(n,2)-toyz(n)*v(n,4))/v(n,1)-qcy(n)
             fez=((v(n,5)+ps(n)-pinfl-tozz(n))*v(n,4) &
                  -toxz(n)*v(n,2)-toyz(n)*v(n,3))/v(n,1)-qcz(n)
!
             sk0= v(n,2)*sn(m,kdir,1) &
                 +v(n,3)*sn(m,kdir,2) &
                 +v(n,4)*sn(m,kdir,3)
             sk1= fxx*sn(m,kdir,1) &
                 +fxy*sn(m,kdir,2) &
                 +fxz*sn(m,kdir,3)
             sk2= fxy*sn(m,kdir,1) &
                 +fyy*sn(m,kdir,2) &
                 +fyz*sn(m,kdir,3)
             sk3= fxz*sn(m,kdir,1) &
                 +fyz*sn(m,kdir,2) &
                 +fzz*sn(m,kdir,3)
             sk4= fex*sn(m,kdir,1) &
                 +fey*sn(m,kdir,2) &
                 +fez*sn(m,kdir,3)
             u(n-ninc,1)=u(n-ninc,1)+sk0
             u(n-ninc,2)=u(n-ninc,2)+sk1
             u(n-ninc,3)=u(n-ninc,3)+sk2
             u(n-ninc,4)=u(n-ninc,4)+sk3
             u(n-ninc,5)=u(n-ninc,5)+sk4
          enddo
       enddo
!
    endif
!
!-----calcul de la 'forcing function'---------------------------
!
    if(ityprk.ne.0) then
       do k=k1,k2m1
          do j=j1,j2m1
             ind1=indc(i1,j,k)
             ind2=indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                ff(n,1) = ff(n,1) - u(n,1)
                ff(n,2) = ff(n,2) - u(n,2)
                ff(n,3) = ff(n,3) - u(n,3)
                ff(n,4) = ff(n,4) - u(n,4)
                ff(n,5) = ff(n,5) - u(n,5)
             enddo
          enddo
       enddo
    endif

    DEALLOCATE(r1,r2,r3,r4,r5)

    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
    end function indc
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
    function    fmp(aa)
      implicit none
      double precision ::  aa,fmp
      fmp=0.25*(1.+sign(1.D0,abs(aa)-1.))*(aa+abs(aa)) &
        +0.125*(1.-sign(1.D0,abs(aa)-1.))*(aa+1.)**2 !*(1.+0.5*(aa-1.)**2)
    end function fmp
    function    fmm(xa)
      implicit none
      double precision :: fmm, xa
      fmm=0.25*(1.+sign(1.D0,abs(xa)-1.))*(xa-abs(xa)) &
        -0.125*(1.-sign(1.D0,abs(xa)-1.))*(xa-1.)**2 !*(1.+0.5*(xa+1.)**2)
    end function fmm
    function    fpp(ta)
      implicit none
      double precision :: fpp, ta
      fpp=0.25*(1.+sign(1.D0,abs(ta)-1.))*(1.+sign(1.D0,abs(ta))) &
        +0.125*(1.-sign(1.D0,abs(ta)-1.))*(ta+1.)**2*(2.-ta+0.75*ta*(ta-1.)**2)
    end function fpp
    function    fpm(ra)
      implicit none
      double precision :: fpm, ra
      fpm=0.25*(1.+sign(1.D0,abs(ra)-1.))*(1.-sign(1.D0,abs(ra))) &
        +0.125*(1.-sign(1.D0,abs(ra)-1.))*(ra-1.)**2*(2.+ra-0.75*ra*(ra+1.)**2)
    end function fpm
    function    phi(a)
      implicit none
      double precision ::   a,phi
      phi=sign(1.D0,a)*max(0.,min(abs(a),sign(1.D0,a)))
    end function phi

  end subroutine sch_ausmp_pond
end module mod_sch_ausmp_pond
