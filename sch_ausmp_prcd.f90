      subroutine sch_ausmp_prcd( &
                 lm,ityprk, &
                 u,v,d,ff, &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 equat, &
                 sn,lgsnlt, &
                 rhol,ul,vl,wl,pl,rhor,ur,vr,wr,prr, &
                 ps)
!
!***********************************************************************
!
!_DA  DATE_C : decembre 2009 - Eric Goncalves / LEGI
!
!     ACT
!_A    Calcul des bilans de flux physiques
!_A    Schema AUSM+ preconditionnee
!_A    Ordre 1, 2 et 3 - extrapolation MUSCL.
!
!***********************************************************************
!
      use para_var
      use para_fige
      use maillage
      use schemanum
      use proprieteflu
      use definition
implicit none
double precision :: fmp
double precision :: fpm
double precision :: fpp
double precision :: fmm
integer :: inc
integer :: indc
double precision :: xa
double precision :: aa
double precision :: ra
double precision :: ta
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: lm
integer :: ityprk
double precision :: u
double precision :: v
double precision :: d
double precision :: ff
double precision :: toxx
double precision :: toxy
double precision :: toxz
double precision :: toyy
double precision :: toyz
double precision :: tozz
double precision :: qcx
double precision :: qcy
double precision :: qcz
double precision :: sn
integer :: lgsnlt
double precision :: rhol
double precision :: ul
double precision :: vl
double precision :: wl
double precision :: pl
double precision :: rhor
double precision :: ur
double precision :: vr
double precision :: wr
double precision :: prr
double precision :: ps
double precision :: fex
double precision :: fey
double precision :: fez
double precision :: fxx
double precision :: fxy
double precision :: fxz
double precision :: fyy
double precision :: fyz
double precision :: fzz
double precision :: hk1
double precision :: hk2
double precision :: hk3
double precision :: hk4
double precision :: hk5
integer :: i1
integer :: i1m1
integer :: i1p1
integer :: i2
integer :: i2m1
integer :: ind1
integer :: ind2
integer :: j1
integer :: j1m1
integer :: j1p1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k1m1
integer :: k1p1
integer :: k2
integer :: k2m1
integer :: kdir
integer :: m
integer :: n
integer :: n0c
integer :: n1
integer :: nci
integer :: ncj
integer :: nck
integer :: nid
integer :: nijd
integer :: ninc
integer :: njd
double precision :: si1
double precision :: si2
double precision :: si3
double precision :: si4
double precision :: si5
double precision :: sj1
double precision :: sj2
double precision :: sj3
double precision :: sj4
double precision :: sj5
double precision :: sk1
double precision :: sk2
double precision :: sk3
double precision :: sk4
double precision :: sk5
!          
!-----------------------------------------------------------------------
!
      character(len=7 ) :: equat
      real nx,ny,nz,cnds,cmm,fpr,dm,qinf,ca,cb,ai,ani,cmi,psi
      real al,ar,q2l,q2r,hl,hr,vnl,vnr,cml,cmr,beta2l,beta2r,betai
      real hi1,hi2,hi3,hi4,hi5,hj1,hj2,hj3,hj4,hj5
      dimension v(ip11,ip60),d(ip11,ip60),u(ip11,ip60),ff(ip11,ip60)
      dimension toxx(ip12),toxy(ip12),toxz(ip12),toyy(ip12),toyz(ip12), &
                tozz(ip12),qcx(ip12),qcy(ip12),qcz(ip12)
      dimension sn(lgsnlt,nind,ndir)
      dimension ps(ip11)
      dimension rhol(ip00),ul(ip00),vl(ip00),wl(ip00),pl(ip00), &
             rhor(ip00),ur(ip00),vr(ip00),wr(ip00),prr(ip00)
!
      indc(i,j,k)=n0c+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
!     fonction du schema ausm+
      fmp(aa)=0.25*(1.+sign(1.,abs(aa)-1.))*(aa+abs(aa)) &
              +0.5*(1.-sign(1.,abs(aa)-1.))* &
             (0.25*(aa+1.)**2+cb*(aa**2-1.)**2)
      fmm(xa)=0.25*(1.+sign(1.,abs(xa)-1.))*(xa-abs(xa)) &
              -0.5*(1.-sign(1.,abs(xa)-1.))* &
             (0.25*(xa-1.)**2-cb**(xa**2-1.)**2)
      fpp(ta)=0.25*(1.+sign(1.,abs(ta)-1.))*(1.+sign(1.,abs(ta))) &
              +0.5*(1.-sign(1.,abs(ta)-1.))* &
             (0.25*(ta+1.)**2*(2.-ta)+ca*ta*(ta**2-1.)**2)
      fpm(ra)=0.25*(1.+sign(1.,abs(ra)-1.))*(1.-sign(1.,abs(ra))) &
              +0.5*(1.-sign(1.,abs(ra)-1.))* &
             (0.25*(ra-1.)**2*(2.+ra)-ca*ra*(ra**2-1.)**2)

      n0c=npc(lm)
      i1=ii1(lm)
      i2=ii2(lm)
      j1=jj1(lm)
      j2=jj2(lm)
      k1=kk1(lm)
      k2=kk2(lm)
!
      nid = id2(lm)-id1(lm)+1
      njd = jd2(lm)-jd1(lm)+1
      nijd = nid*njd
!
      i1p1=i1+1
      j1p1=j1+1
      k1p1=k1+1
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
!
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
      qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
      ca=0.
      cb=0.
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
!-----definition des variables extrapolees--------------------------------
!
      do k=k1,k2m1
       do j=j1,j2m1
        ind1 = indc(i1p1,j,k)
        ind2 = indc(i2m1,j,k)
        do n=ind1,ind2
         m=n-n0c
         rhol(m)=v(n-ninc,1)+0.25*muscl*( &
           (1.-xk)*(v(n-ninc,1)-v(n-2*ninc,1)) &
          +(1.+xk)*(v(n     ,1)-v(n-ninc  ,1)))
         ul(m)=v(n-ninc,2)/v(n-ninc,1)+0.25*muscl*( &
          (1.-xk)*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
         +(1.+xk)*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
         vl(m)=v(n-ninc,3)/v(n-ninc,1)+0.25*muscl*( &
          (1.-xk)*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
         +(1.+xk)*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
         wl(m)=v(n-ninc,4)/v(n-ninc,1)+0.25*muscl*( &
          (1.-xk)*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
         +(1.+xk)*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
         pl(m)=ps(n-ninc)+0.25*muscl*( &
           (1.-xk)*(ps(n-ninc)-ps(n-2*ninc)) &
          +(1.+xk)*(ps(n     )-ps(n-  ninc)))
!
         rhor(m)=v(n,1)-0.25*muscl*((1.+xk)*(v(n,1)     -v(n-ninc,1)) &
                                   +(1.-xk)*(v(n+ninc,1)-v(n     ,1)))
         ur(m)=v(n,2)/v(n,1)-0.25*muscl*( &
           (1.+xk)*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
          +(1.-xk)*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
         vr(m)=v(n,3)/v(n,1)-0.25*muscl*( &
           (1.+xk)*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
          +(1.-xk)*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
         wr(m)=v(n,4)/v(n,1)-0.25*muscl*( &
           (1.+xk)*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
          +(1.-xk)*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
         prr(m)=ps(n)-0.25*muscl*((1.+xk)*(ps(n)     -ps(n-ninc)) &
                                +(1.-xk)*(ps(n+ninc)-ps(n     )))
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
         cnds=sqrt(sn(m,kdir,1)*sn(m,kdir,1)+ &
                   sn(m,kdir,2)*sn(m,kdir,2)+ &
                   sn(m,kdir,3)*sn(m,kdir,3))
         nx=sn(m,kdir,1)/cnds
         ny=sn(m,kdir,2)/cnds
         nz=sn(m,kdir,3)/cnds
!        calcul des etats gauche et droit
         al=sqrt(gam*pl(m)/rhol(m))
         ar=sqrt(gam*prr(m)/rhor(m))
         q2l=ul(m)**2+vl(m)**2+wl(m)**2
         q2r=ur(m)**2+vr(m)**2+wr(m)**2
         hl=al*al/gam1+0.5*q2l
         hr=ar*ar/gam1+0.5*q2r
         vnl=ul(m)*nx+vl(m)*ny+wl(m)*nz
         vnr=ur(m)*nx+vr(m)*ny+wr(m)*nz
!        calcul de la vitesse du son a l'interface
         ai=sqrt(al*ar)
!        terme de preconditionnement
         beta2l=min(max(q2l/al**2,cte*(qinf/al)**2),1.)
         beta2r=min(max(q2r/ar**2,cte*(qinf/ar)**2),1.)
         betai=0.5*(sqrt(beta2l)+sqrt(beta2r))
         cmm=0.5*(vnl+vnr)/ai
         fpr=sqrt((1.-betai**2)**2*cmm**2+4.*betai**2)/(1.+betai**2)
!        calcul de la vitesse du son numerique a l'interface
         ani=fpr*ai
!        calcul des nombres de Mach gauche et droit
         cml=0.5*((1.+betai**2)*vnl+(1.-betai**2)*vnr)/ani
         cmr=0.5*((1.+betai**2)*vnr+(1.-betai**2)*vnl)/ani
!        calcul du nombre de Mach a l'interface
         cmi=fmp(cml)+fmm(cmr)
!        calcul de la pression statique a l'interface (x2)
         psi=2.*(pl(m)*fpp(cml)+prr(m)*fpm(cmr))-2.*pinfl
!        calcul du flux de masse
         dm=0.5*ani*(cmi*(rhol(m)+rhor(m))-abs(cmi)*(rhor(m)-rhol(m)))*cnds
!        calcul du flux numerique a l'interface
         hi1=dm
         hi2=dm*(ul(m)+ur(m))-abs(dm)*(ur(m)-ul(m))+psi*sn(m,kdir,1) &
            - (toxx(n)+toxx(n-ninc))*sn(m,kdir,1) &
            - (toxy(n)+toxy(n-ninc))*sn(m,kdir,2) &
            - (toxz(n)+toxz(n-ninc))*sn(m,kdir,3)
         hi3=dm*(vl(m)+vr(m))-abs(dm)*(vr(m)-vl(m))+psi*sn(m,kdir,2) &
            - (toxy(n)+toxy(n-ninc))*sn(m,kdir,1) &
            - (toyy(n)+toyy(n-ninc))*sn(m,kdir,2) &
            - (toyz(n)+toyz(n-ninc))*sn(m,kdir,3)
         hi4=dm*(wl(m)+wr(m))-abs(dm)*(wr(m)-wl(m))+psi*sn(m,kdir,3) &
            - (toxz(n)+toxz(n-ninc))*sn(m,kdir,1) &
            - (toyz(n)+toyz(n-ninc))*sn(m,kdir,2) &
            - (tozz(n)+tozz(n-ninc))*sn(m,kdir,3)
         hi5=dm*(hl+hr)-abs(dm)*(hr-hl) &
            - (d(n,1)+d(n-ninc,1))*sn(m,kdir,1) &
            - (d(n,2)+d(n-ninc,2))*sn(m,kdir,2) &
            - (d(n,3)+d(n-ninc,3))*sn(m,kdir,3)
!
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
!DEC$ IVDEP
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
        si1= v(n-ninc,2)*sn(m,kdir,1) &
            +v(n-ninc,3)*sn(m,kdir,2) &
            +v(n-ninc,4)*sn(m,kdir,3)
        si2= fxx*sn(m,kdir,1) &
            +fxy*sn(m,kdir,2) &
            +fxz*sn(m,kdir,3)
        si3= fxy*sn(m,kdir,1) &
            +fyy*sn(m,kdir,2) &
            +fyz*sn(m,kdir,3)
        si4= fxz*sn(m,kdir,1) &
            +fyz*sn(m,kdir,2) &
            +fzz*sn(m,kdir,3)
        si5= fex*sn(m,kdir,1) &
            +fey*sn(m,kdir,2) &
            +fez*sn(m,kdir,3)
        u(n,1)=u(n,1)-si1
        u(n,2)=u(n,2)-si2
        u(n,3)=u(n,3)-si3
        u(n,4)=u(n,4)-si4
        u(n,5)=u(n,5)-si5
       enddo
      enddo
!
      do k=k1,k2m1
       ind1 = indc(i2,j1  ,k)
       ind2 = indc(i2,j2m1,k)
!DEC$ IVDEP
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
        si1= v(n,2)*sn(m,kdir,1) &
            +v(n,3)*sn(m,kdir,2) &
            +v(n,4)*sn(m,kdir,3)
        si2= fxx*sn(m,kdir,1) &
            +fxy*sn(m,kdir,2) &
            +fxz*sn(m,kdir,3)
        si3= fxy*sn(m,kdir,1) &
            +fyy*sn(m,kdir,2) &
            +fyz*sn(m,kdir,3)
        si4= fxz*sn(m,kdir,1) &
            +fyz*sn(m,kdir,2) &
            +fzz*sn(m,kdir,3)
        si5= fex*sn(m,kdir,1) &
            +fey*sn(m,kdir,2) &
            +fez*sn(m,kdir,3)
        u(n-ninc,1)=u(n-ninc,1)+si1
        u(n-ninc,2)=u(n-ninc,2)+si2
        u(n-ninc,3)=u(n-ninc,3)+si3
        u(n-ninc,4)=u(n-ninc,4)+si4
        u(n-ninc,5)=u(n-ninc,5)+si5
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
      do k=k1,k2m1
       do j=j1p1,j2m1
        ind1 = indc(i1  ,j,k)
        ind2 = indc(i2m1,j,k)
        do n=ind1,ind2
         m=n-n0c
         rhol(m)=v(n-ninc,1)+0.25*muscl*( &
           (1.-xk)*(v(n-ninc,1)-v(n-2*ninc,1)) &
          +(1.+xk)*(v(n     ,1)-v(n-ninc  ,1)))
         ul(m)=v(n-ninc,2)/v(n-ninc,1)+0.25*muscl*( &
          (1.-xk)*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
         +(1.+xk)*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
         vl(m)=v(n-ninc,3)/v(n-ninc,1)+0.25*muscl*( &
          (1.-xk)*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
         +(1.+xk)*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
         wl(m)=v(n-ninc,4)/v(n-ninc,1)+0.25*muscl*( &
          (1.-xk)*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
         +(1.+xk)*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
         pl(m)=ps(n-ninc)+0.25*muscl*( &
           (1.-xk)*(ps(n-ninc)-ps(n-2*ninc)) &
          +(1.+xk)*(ps(n     )-ps(n-  ninc)))
!
         rhor(m)=v(n,1)-0.25*muscl*((1.+xk)*(v(n,1)     -v(n-ninc,1)) &
                                   +(1.-xk)*(v(n+ninc,1)-v(n     ,1)))
         ur(m)=v(n,2)/v(n,1)-0.25*muscl*( &
           (1.+xk)*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
          +(1.-xk)*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
         vr(m)=v(n,3)/v(n,1)-0.25*muscl*( &
           (1.+xk)*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
          +(1.-xk)*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
         wr(m)=v(n,4)/v(n,1)-0.25*muscl*( &
           (1.+xk)*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
          +(1.-xk)*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
         prr(m)=ps(n)-0.25*muscl*((1.+xk)*(ps(n)     -ps(n-ninc)) &
                                +(1.-xk)*(ps(n+ninc)-ps(n     )))
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
         cnds=sqrt(sn(m,kdir,1)*sn(m,kdir,1)+ &
                   sn(m,kdir,2)*sn(m,kdir,2)+ &
                   sn(m,kdir,3)*sn(m,kdir,3))
         nx=sn(m,kdir,1)/cnds
         ny=sn(m,kdir,2)/cnds
         nz=sn(m,kdir,3)/cnds
!        calcul des etats gauche et droit
         al=sqrt(gam*pl(m)/rhol(m))
         ar=sqrt(gam*prr(m)/rhor(m))
         q2l=ul(m)**2+vl(m)**2+wl(m)**2
         q2r=ur(m)**2+vr(m)**2+wr(m)**2
         hl=al*al/gam1+0.5*q2l
         hr=ar*ar/gam1+0.5*q2r
         vnl=ul(m)*nx+vl(m)*ny+wl(m)*nz
         vnr=ur(m)*nx+vr(m)*ny+wr(m)*nz
!        calcul de la vitesse du son a l'interface
         ai=sqrt(al*ar)
!        terme de preconditionnement
         beta2l=min(max(q2l/al**2,cte*(qinf/al)**2),1.)
         beta2r=min(max(q2r/ar**2,cte*(qinf/ar)**2),1.)
         betai=0.5*(sqrt(beta2l)+sqrt(beta2r))
         cmm=0.5*(vnl+vnr)/ai
         fpr=sqrt((1.-betai**2)**2*cmm**2+4.*betai**2)/(1.+betai**2)
!        calcul de la vitesse du son numerique a l'interface
         ani=fpr*ai
!        calcul des nombres de Mach gauche et droit
         cml=0.5*((1.+betai**2)*vnl+(1.-betai**2)*vnr)/ani
         cmr=0.5*((1.+betai**2)*vnr+(1.-betai**2)*vnl)/ani
!        calcul du nombre de Mach a l'interface
         cmi=fmp(cml)+fmm(cmr)
!        calcul de la pression statique a l'interface (x2)
         psi=2.*(pl(m)*fpp(cml)+prr(m)*fpm(cmr))-2.*pinfl
!        calcul du flux de masse
         dm=0.5*ani*(cmi*(rhol(m)+rhor(m))-abs(cmi)*(rhor(m)-rhol(m)))*cnds
!        calcul du flux numerique
         hj1=dm
         hj2=dm*(ul(m)+ur(m))-abs(dm)*(ur(m)-ul(m))+psi*sn(m,kdir,1) &
            - (toxx(n)+toxx(n-ninc))*sn(m,kdir,1) &
            - (toxy(n)+toxy(n-ninc))*sn(m,kdir,2) &
            - (toxz(n)+toxz(n-ninc))*sn(m,kdir,3)
         hj3=dm*(vl(m)+vr(m))-abs(dm)*(vr(m)-vl(m))+psi*sn(m,kdir,2) &
            - (toxy(n)+toxy(n-ninc))*sn(m,kdir,1) &
            - (toyy(n)+toyy(n-ninc))*sn(m,kdir,2) &
            - (toyz(n)+toyz(n-ninc))*sn(m,kdir,3)
         hj4=dm*(wl(m)+wr(m))-abs(dm)*(wr(m)-wl(m))+psi*sn(m,kdir,3) &
            - (toxz(n)+toxz(n-ninc))*sn(m,kdir,1) &
            - (toyz(n)+toyz(n-ninc))*sn(m,kdir,2) &
            - (tozz(n)+tozz(n-ninc))*sn(m,kdir,3)
         hj5=dm*(hl+hr)-abs(dm)*(hr-hl) &
            - (d(n,1)+d(n-ninc,1))*sn(m,kdir,1) &
            - (d(n,2)+d(n-ninc,2))*sn(m,kdir,2) &
            - (d(n,3)+d(n-ninc,3))*sn(m,kdir,3)
!
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
!DEC$ IVDEP
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
        sj1= v(n-ninc,2)*sn(m,kdir,1) &
            +v(n-ninc,3)*sn(m,kdir,2) &
            +v(n-ninc,4)*sn(m,kdir,3)
        sj2= fxx*sn(m,kdir,1) &
            +fxy*sn(m,kdir,2) &
            +fxz*sn(m,kdir,3)
        sj3= fxy*sn(m,kdir,1) &
            +fyy*sn(m,kdir,2) &
            +fyz*sn(m,kdir,3)
        sj4= fxz*sn(m,kdir,1) &
            +fyz*sn(m,kdir,2) &
            +fzz*sn(m,kdir,3)
        sj5= fex*sn(m,kdir,1) &
            +fey*sn(m,kdir,2) &
            +fez*sn(m,kdir,3)
        u(n,1)=u(n,1)-sj1
        u(n,2)=u(n,2)-sj2
        u(n,3)=u(n,3)-sj3
        u(n,4)=u(n,4)-sj4
        u(n,5)=u(n,5)-sj5
       enddo
      enddo
!
      do k=k1,k2m1
       ind1 = indc(i1  ,j2,k)
       ind2 = indc(i2m1,j2,k)
!DEC$ IVDEP
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
        sj1= v(n,2)*sn(m,kdir,1) &
            +v(n,3)*sn(m,kdir,2) &
            +v(n,4)*sn(m,kdir,3)
        sj2= fxx*sn(m,kdir,1) &
            +fxy*sn(m,kdir,2) &
            +fxz*sn(m,kdir,3)
        sj3= fxy*sn(m,kdir,1) &
            +fyy*sn(m,kdir,2) &
            +fyz*sn(m,kdir,3)
        sj4= fxz*sn(m,kdir,1) &
            +fyz*sn(m,kdir,2) &
            +fzz*sn(m,kdir,3)
        sj5= fex*sn(m,kdir,1) &
            +fey*sn(m,kdir,2) &
            +fez*sn(m,kdir,3)
        u(n-ninc,1)=u(n-ninc,1)+sj1
        u(n-ninc,2)=u(n-ninc,2)+sj2
        u(n-ninc,3)=u(n-ninc,3)+sj3
        u(n-ninc,4)=u(n-ninc,4)+sj4
        u(n-ninc,5)=u(n-ninc,5)+sj5
       enddo
      enddo
!
!------direction k----------------------------------------------
!
      if(equat(3:4).eq.'3d') then
       kdir=3
       ninc=nck
!
      do k=k1p1,k2m1
       do j=j1,j2m1
        ind1 = indc(i1,j,k)
        ind2 = indc(i2m1,j,k)
        do n=ind1,ind2
         m=n-n0c
         rhol(m)=v(n-ninc,1)+0.25*muscl*( &
           (1.-xk)*(v(n-ninc,1)-v(n-2*ninc,1)) &
          +(1.+xk)*(v(n     ,1)-v(n-ninc  ,1)))
         ul(m)=v(n-ninc,2)/v(n-ninc,1)+0.25*muscl*( &
          (1.-xk)*(v(n-ninc,2)/v(n-ninc,1)-v(n-2*ninc,2)/v(n-2*ninc,1)) &
         +(1.+xk)*(v(n     ,2)/v(n     ,1)-v(n-ninc  ,2)/v(n-ninc  ,1)))
         vl(m)=v(n-ninc,3)/v(n-ninc,1)+0.25*muscl*( &
          (1.-xk)*(v(n-ninc,3)/v(n-ninc,1)-v(n-2*ninc,3)/v(n-2*ninc,1)) &
         +(1.+xk)*(v(n     ,3)/v(n     ,1)-v(n-ninc  ,3)/v(n-ninc  ,1)))
         wl(m)=v(n-ninc,4)/v(n-ninc,1)+0.25*muscl*( &
          (1.-xk)*(v(n-ninc,4)/v(n-ninc,1)-v(n-2*ninc,4)/v(n-2*ninc,1)) &
         +(1.+xk)*(v(n     ,4)/v(n     ,1)-v(n-ninc  ,4)/v(n-ninc  ,1)))
         pl(m)=ps(n-ninc)+0.25*muscl*( &
           (1.-xk)*(ps(n-ninc)-ps(n-2*ninc)) &
          +(1.+xk)*(ps(n     )-ps(n-  ninc)))
!
         rhor(m)=v(n,1)-0.25*muscl*((1.+xk)*(v(n,1)     -v(n-ninc,1)) &
                                 +(1.-xk)*(v(n+ninc,1)-v(n     ,1)))
         ur(m)=v(n,2)/v(n,1)-0.25*muscl*( &
           (1.+xk)*(v(n     ,2)/v(n     ,1)-v(n-ninc,2)/v(n-ninc,1)) &
          +(1.-xk)*(v(n+ninc,2)/v(n+ninc,1)-v(n     ,2)/v(n     ,1)))
         vr(m)=v(n,3)/v(n,1)-0.25*muscl*( &
           (1.+xk)*(v(n     ,3)/v(n     ,1)-v(n-ninc,3)/v(n-ninc,1)) &
          +(1.-xk)*(v(n+ninc,3)/v(n+ninc,1)-v(n     ,3)/v(n     ,1)))
         wr(m)=v(n,4)/v(n,1)-0.25*muscl*( &
           (1.+xk)*(v(n     ,4)/v(n     ,1)-v(n-ninc,4)/v(n-ninc,1)) &
          +(1.-xk)*(v(n+ninc,4)/v(n+ninc,1)-v(n     ,4)/v(n     ,1)))
         prr(m)=ps(n)-0.25*muscl*((1.+xk)*(ps(n)     -ps(n-ninc)) &
                              +(1.-xk)*(ps(n+ninc)-ps(n     )))
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
         cnds=sqrt(sn(m,kdir,1)*sn(m,kdir,1)+ &
                   sn(m,kdir,2)*sn(m,kdir,2)+ &
                   sn(m,kdir,3)*sn(m,kdir,3))
         nx=sn(m,kdir,1)/cnds
         ny=sn(m,kdir,2)/cnds
         nz=sn(m,kdir,3)/cnds
!        calcul des etats gauche et droit
         al=sqrt(gam*pl(m)/rhol(m))
         ar=sqrt(gam*prr(m)/rhor(m))
         q2l=ul(m)**2+vl(m)**2+wl(m)**2
         q2r=ur(m)**2+vr(m)**2+wr(m)**2
         hl=al*al/gam1+0.5*q2l
         hr=ar*ar/gam1+0.5*q2r
         vnl=ul(m)*nx+vl(m)*ny+wl(m)*nz
         vnr=ur(m)*nx+vr(m)*ny+wr(m)*nz
!        calcul de la vitesse du son a l'interface
         ai=sqrt(al*ar)
!        terme de preconditionnement
         beta2l=min(max(q2l/al**2,cte*(qinf/al)**2),1.)
         beta2r=min(max(q2r/ar**2,cte*(qinf/ar)**2),1.)
         betai=0.5*(sqrt(beta2l)+sqrt(beta2r))
         cmm=0.5*(vnl+vnr)/ai
         fpr=sqrt((1.-betai**2)**2*cmm**2+4.*betai**2)/(1.+betai**2)
!        calcul de la vitesse du son numerique a l'interface
         ani=fpr*ai
!        calcul des nombres de Mach gauche et droit
         cml=0.5*((1.+betai**2)*vnl+(1.-betai**2)*vnr)/ani
         cmr=0.5*((1.+betai**2)*vnr+(1.-betai**2)*vnl)/ani
!        calcul du nombre de Mach a l'interface
         cmi=fmp(cml)+fmm(cmr)
!        calcul de la pression statique a l'interface (x2)
         psi=2.*(pl(m)*fpp(cml)+prr(m)*fpm(cmr))-2.*pinfl
!        calcul du flux de masse
         dm=0.5*ani*(cmi*(rhol(m)+rhor(m))-abs(cmi)*(rhor(m)-rhol(m)))*cnds
!        calcul du flux numerique
         hk1=dm
         hk2=dm*(ul(m)+ur(m))-abs(dm)*(ur(m)-ul(m))+psi*sn(m,kdir,1) &
            - (toxx(n)+toxx(n-ninc))*sn(m,kdir,1) &
            - (toxy(n)+toxy(n-ninc))*sn(m,kdir,2) &
            - (toxz(n)+toxz(n-ninc))*sn(m,kdir,3)
         hk3=dm*(vl(m)+vr(m))-abs(dm)*(vr(m)-vl(m))+psi*sn(m,kdir,2) &
            - (toxy(n)+toxy(n-ninc))*sn(m,kdir,1) &
            - (toyy(n)+toyy(n-ninc))*sn(m,kdir,2) &
            - (toyz(n)+toyz(n-ninc))*sn(m,kdir,3)
         hk4=dm*(wl(m)+wr(m))-abs(dm)*(wr(m)-wl(m))+psi*sn(m,kdir,3) &
            - (toxz(n)+toxz(n-ninc))*sn(m,kdir,1) &
            - (toyz(n)+toyz(n-ninc))*sn(m,kdir,2) &
            - (tozz(n)+tozz(n-ninc))*sn(m,kdir,3)
         hk5=dm*(hl+hr)-abs(dm)*(hr-hl) &
            - (d(n,1)+d(n-ninc,1))*sn(m,kdir,1) &
            - (d(n,2)+d(n-ninc,2))*sn(m,kdir,2) &
            - (d(n,3)+d(n-ninc,3))*sn(m,kdir,3)
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
!DEC$ IVDEP
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
        sk1= v(n-ninc,2)*sn(m,kdir,1) &
            +v(n-ninc,3)*sn(m,kdir,2) &
            +v(n-ninc,4)*sn(m,kdir,3)
        sk2= fxx*sn(m,kdir,1) &
            +fxy*sn(m,kdir,2) &
            +fxz*sn(m,kdir,3)
        sk3= fxy*sn(m,kdir,1) &
            +fyy*sn(m,kdir,2) &
            +fyz*sn(m,kdir,3)
        sk4= fxz*sn(m,kdir,1) &
            +fyz*sn(m,kdir,2) &
            +fzz*sn(m,kdir,3)
        sk5= fex*sn(m,kdir,1) &
            +fey*sn(m,kdir,2) &
            +fez*sn(m,kdir,3)
        u(n,1)=u(n,1)-sk1
        u(n,2)=u(n,2)-sk2
        u(n,3)=u(n,3)-sk3
        u(n,4)=u(n,4)-sk4
        u(n,5)=u(n,5)-sk5
       enddo
      enddo
!
      do j=j1,j2m1
       ind1 = indc(i1  ,j,k2)
       ind2 = indc(i2m1,j,k2)
!DEC$ IVDEP
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
        sk1= v(n,2)*sn(m,kdir,1) &
            +v(n,3)*sn(m,kdir,2) &
            +v(n,4)*sn(m,kdir,3)
        sk2= fxx*sn(m,kdir,1) &
            +fxy*sn(m,kdir,2) &
            +fxz*sn(m,kdir,3)
        sk3= fxy*sn(m,kdir,1) &
            +fyy*sn(m,kdir,2) &
            +fyz*sn(m,kdir,3)
        sk4= fxz*sn(m,kdir,1) &
            +fyz*sn(m,kdir,2) &
            +fzz*sn(m,kdir,3)
        sk5= fex*sn(m,kdir,1) &
            +fey*sn(m,kdir,2) &
            +fez*sn(m,kdir,3)
        u(n-ninc,1)=u(n-ninc,1)+sk1
        u(n-ninc,2)=u(n-ninc,2)+sk2
        u(n-ninc,3)=u(n-ninc,3)+sk3
        u(n-ninc,4)=u(n-ninc,4)+sk4
        u(n-ninc,5)=u(n-ninc,5)+sk5
       enddo
      enddo
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

      return
      end
