module mod_sch_rusanov_prcd
implicit none
contains
      subroutine sch_rusanov_prcd( &
                 lm,ityprk, &
                 u,v,ff, &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 equat, &
                 sn,lgsnlt, &
                 rhol,ul,vl,wl,pl,rhor,ur,vr,wr,prr, &
                 ps)
!
!***********************************************************************
!
!_DA  DATE_C : octobre 2011 - Eric Goncalves / LEGI
!
!     ACT
!_A    Calcul des bilans de flux physiques. Schema de Rusanov.
!_A    Avec extrapolation MUSCL (pour ordre 2 et 3).
!_A    Pas de limiteur de pentes
!_A    Preconditionnement basse vitesse de Turkel (P,u,e).
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use proprieteflu
      use schemanum
      use definition
implicit none
integer :: inc
integer :: indc
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
double precision :: al
double precision :: ar
double precision :: beta2l
double precision :: beta2r
double precision :: cnds
double precision :: dfex
double precision :: dfey
double precision :: dfez
double precision :: dfxx
double precision :: dfxy
double precision :: dfxz
double precision :: dfyy
double precision :: dfyz
double precision :: dfzz
double precision :: el
double precision :: er
double precision :: fex
double precision :: fey
double precision :: fez
double precision :: fxx
double precision :: fxy
double precision :: fxz
double precision :: fyy
double precision :: fyz
double precision :: fzz
double precision :: hl
double precision :: hr
integer :: i1
integer :: i1m1
integer :: i1p1
integer :: i2
integer :: i2m1
integer :: ind1
integer :: ind2
integer :: isortie
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
double precision :: q2l
double precision :: q2r
double precision :: qinf
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
double precision :: vnl
double precision :: vnr
!
!-----------------------------------------------------------------------
!
      real nx,ny,nz,ids
      real f1,f2,f3,f4,f5,g1,g2,g3,g4,g5,h1,h2,h3,h4,h5,rm,cpl,cpr
      real di1,di2,di3,di4,di5,dj1,dj2,dj3,dj4,dj5,dk1,dk2,dk3,dk4,dk5
      character(len=7 ) :: equat
      dimension v(ip11,ip60),u(ip11,ip60),ff(ip11,ip60)
      dimension toxx(ip12),toxy(ip12),toxz(ip12),toyy(ip12),toyz(ip12), &
                tozz(ip12),qcx(ip12),qcy(ip12),qcz(ip12)
      dimension sn(lgsnlt,nind,ndir)
      dimension ps(ip11)
      dimension rhol(ip00),ul(ip00),vl(ip00),wl(ip00),pl(ip00), &
             rhor(ip00),ur(ip00),vr(ip00),wr(ip00),prr(ip00)
!
      indc(i,j,k)=n0c+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd

      isortie=0
!
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
      qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
!
!-----calcul des densites de flux convectives -----------------------------
!
      if(equat(3:5).eq.'2dk') then
       ind1 = indc(i1m1,j1m1,k1  )
       ind2 = indc(i2  ,j2  ,k2m1)
      elseif(equat(3:4).eq.'3d') then
       ind1 = indc(i1m1,j1m1,k1m1)
       ind2 = indc(i2  ,j2  ,k2  )
      endif
       do n=ind1,ind2
        u(n,1)=0.
        u(n,2)=0.
        u(n,3)=0.
        u(n,4)=0.
        u(n,5)=0.
       enddo
!
!*******************************************************************************
! calcul du flux numerique par direction suivant les etapes successives :
!    1) evaluation des variables primitives extrapolees
!    2) evaluation du flux numerique
!*******************************************************************************
!
!------direction i-------------------------------------------------------
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
!DEC$ IVDEP
        do n=ind1,ind2
         m=n-n0c
!        vecteur normal unitaire a la face consideree
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
         el=pl(m)/(gam1*rhol(m))+0.5*q2l+pinfl/rhol(m)
         er=prr(m)/(gam1*rhor(m))+0.5*q2r+pinfl/rhor(m)
         vnl=ul(m)*nx+vl(m)*ny+wl(m)*nz
         vnr=ur(m)*nx+vr(m)*ny+wr(m)*nz
!        calcul des vitesses caracteristiques preconditionnees
         beta2l=min(max(q2l/al**2,cte*(qinf/al)**2),1.)
         beta2r=min(max(q2r/ar**2,cte*(qinf/ar)**2),1.)
         cpl=0.5*((1.+beta2l)*vnl*cnds+sqrt(((1.-beta2l)*vnl)**2+4.*beta2l*cnds**2*al**2))
         cpr=0.5*((1.+beta2r)*vnr*cnds+sqrt(((1.-beta2r)*vnr)**2+4.*beta2r*cnds**2*ar**2))
!        calcul rayon spectral
         rm=max(abs(cpl),abs(cpr))
!        evaluation du terme de dissipation (x2)
         di1=rm*(rhor(m)-rhol(m))
         di2=rm*(rhor(m)*ur(m)-rhol(m)*ul(m))
         di3=rm*(rhor(m)*vr(m)-rhol(m)*vl(m))
         di4=rm*(rhor(m)*wr(m)-rhol(m)*wl(m))
         di5=rm*(rhor(m)*er-rhol(m)*el)
!        calcul du flux a l'interface i-1/2 (x2)
         dfxx=rhor(m)*ur(m)**2 + prr(m)-pinfl - toxx(n) &
             +rhol(m)*ul(m)**2 + pl(m)-pinfl - toxx(n-ninc)
         dfxy=rhor(m)*ur(m)*vr(m) - toxy(n) &
             +rhol(m)*ul(m)*vl(m) - toxy(n-ninc)
         dfxz=rhor(m)*ur(m)*wr(m) - toxz(n) &
             +rhol(m)*ul(m)*wl(m) - toxz(n-ninc)
         dfyy=rhor(m)*vr(m)**2 + prr(m)-pinfl - toyy(n) &
             +rhol(m)*vl(m)**2 + pl(m)-pinfl - toyy(n-ninc)
         dfyz=rhor(m)*vr(m)*wr(m) - toyz(n) &
             +rhol(m)*vl(m)*wl(m) - toyz(n-ninc)
         dfzz=rhor(m)*wr(m)**2 + prr(m)-pinfl - tozz(n) &
             +rhol(m)*wl(m)**2 + pl(m)-pinfl - tozz(n-ninc)
         dfex=(rhor(m)*er+prr(m)-pinfl)*ur(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*ul(m) &
             -(toxx(n     )*ur(m)+toxy(n     )*vr(m)+toxz(n     )*wr(m)) &
             -(toxx(n-ninc)*ul(m)+toxy(n-ninc)*vl(m)+toxz(n-ninc)*wl(m)) &
             - qcx(n) - qcx(n-ninc)
         dfey=(rhor(m)*er+prr(m)-pinfl)*vr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*vl(m) &
             -(toxy(n     )*ur(m)+toyy(n     )*vr(m)+toyz(n     )*wr(m)) &
             -(toxy(n-ninc)*ul(m)+toyy(n-ninc)*vl(m)+toyz(n-ninc)*wl(m)) &
             - qcy(n) - qcy(n-ninc)
         dfez=(rhor(m)*er+prr(m)-pinfl)*wr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*wl(m) &
             -(toxz(n     )*ur(m)+toyz(n     )*vr(m)+tozz(n     )*wr(m)) &
             -(toxz(n-ninc)*ul(m)+toyz(n-ninc)*vl(m)+tozz(n-ninc)*wl(m)) &
             - qcz(n) - qcz(n-ninc)
         f1=(rhor(m)*ur(m)+rhol(m)*ul(m))*sn(m,kdir,1) &
           +(rhor(m)*vr(m)+rhol(m)*vl(m))*sn(m,kdir,2) &
           +(rhor(m)*wr(m)+rhol(m)*wl(m))*sn(m,kdir,3)-di1
         f2=dfxx*sn(m,kdir,1)+dfxy*sn(m,kdir,2)+dfxz*sn(m,kdir,3)-di2
         f3=dfxy*sn(m,kdir,1)+dfyy*sn(m,kdir,2)+dfyz*sn(m,kdir,3)-di3
         f4=dfxz*sn(m,kdir,1)+dfyz*sn(m,kdir,2)+dfzz*sn(m,kdir,3)-di4
         f5=dfex*sn(m,kdir,1)+dfey*sn(m,kdir,2)+dfez*sn(m,kdir,3)-di5
         u(n,1)=u(n,1)-0.5*f1
         u(n,2)=u(n,2)-0.5*f2
         u(n,3)=u(n,3)-0.5*f3
         u(n,4)=u(n,4)-0.5*f4
         u(n,5)=u(n,5)-0.5*f5
         u(n-ninc,1)=u(n-ninc,1)+0.5*f1
         u(n-ninc,2)=u(n-ninc,2)+0.5*f2
         u(n-ninc,3)=u(n-ninc,3)+0.5*f3
         u(n-ninc,4)=u(n-ninc,4)+0.5*f4
         u(n-ninc,5)=u(n-ninc,5)+0.5*f5
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
!DEC$ IVDEP
        do n=ind1,ind2
         m=n-n0c
!        vecteur normal unitaire a la face consideree
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
         el=pl(m)/(gam1*rhol(m))+0.5*q2l+pinfl/rhol(m)
         er=prr(m)/(gam1*rhor(m))+0.5*q2r+pinfl/rhor(m)
         vnl=ul(m)*nx+vl(m)*ny+wl(m)*nz
         vnr=ur(m)*nx+vr(m)*ny+wr(m)*nz
!        calcul des vitesses caracteristiques preconditionnees
         beta2l=min(max(q2l/al**2,cte*(qinf/al)**2),1.)
         beta2r=min(max(q2r/ar**2,cte*(qinf/ar)**2),1.)
         cpl=0.5*((1.+beta2l)*vnl*cnds+sqrt(((1.-beta2l)*vnl)**2+4.*beta2l*cnds**2*al**2))
         cpr=0.5*((1.+beta2r)*vnr*cnds+sqrt(((1.-beta2r)*vnr)**2+4.*beta2r*cnds**2*ar**2))
!        calcul rayon spectral
         rm=max(abs(cpl),abs(cpr))
!        evaluation du terme de dissipation (x2)
         dj1=rm*(rhor(m)-rhol(m))
         dj2=rm*(rhor(m)*ur(m)-rhol(m)*ul(m))
         dj3=rm*(rhor(m)*vr(m)-rhol(m)*vl(m))
         dj4=rm*(rhor(m)*wr(m)-rhol(m)*wl(m))
         dj5=rm*(rhor(m)*er-rhol(m)*el)
!        calcul du flux a l'interface j-1/2 (x2)
         dfxx=rhor(m)*ur(m)**2 + prr(m)-pinfl - toxx(n) &
             +rhol(m)*ul(m)**2 + pl(m)-pinfl - toxx(n-ninc)
         dfxy=rhor(m)*ur(m)*vr(m) - toxy(n) &
             +rhol(m)*ul(m)*vl(m) - toxy(n-ninc)
         dfxz=rhor(m)*ur(m)*wr(m) - toxz(n) &
             +rhol(m)*ul(m)*wl(m) - toxz(n-ninc)
         dfyy=rhor(m)*vr(m)**2 + prr(m)-pinfl - toyy(n) &
             +rhol(m)*vl(m)**2 + pl(m)-pinfl - toyy(n-ninc)
         dfyz=rhor(m)*vr(m)*wr(m) - toyz(n) &
             +rhol(m)*vl(m)*wl(m) - toyz(n-ninc)
         dfzz=rhor(m)*wr(m)**2 + prr(m)-pinfl - tozz(n) &
             +rhol(m)*wl(m)**2 + pl(m)-pinfl - tozz(n-ninc)
         dfex=(rhor(m)*er+prr(m)-pinfl)*ur(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*ul(m) &
             -(toxx(n     )*ur(m)+toxy(n     )*vr(m)+toxz(n     )*wr(m)) &
             -(toxx(n-ninc)*ul(m)+toxy(n-ninc)*vl(m)+toxz(n-ninc)*wl(m)) &
             - qcx(n) - qcx(n-ninc)
         dfey=(rhor(m)*er+prr(m)-pinfl)*vr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*vl(m) &
             -(toxy(n     )*ur(m)+toyy(n     )*vr(m)+toyz(n     )*wr(m)) &
             -(toxy(n-ninc)*ul(m)+toyy(n-ninc)*vl(m)+toyz(n-ninc)*wl(m)) &
             - qcy(n) - qcy(n-ninc)
         dfez=(rhor(m)*er+prr(m)-pinfl)*wr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*wl(m) &
             -(toxz(n     )*ur(m)+toyz(n     )*vr(m)+tozz(n     )*wr(m)) &
             -(toxz(n-ninc)*ul(m)+toyz(n-ninc)*vl(m)+tozz(n-ninc)*wl(m)) &
             - qcz(n) - qcz(n-ninc)
!
         g1=(rhor(m)*ur(m)+rhol(m)*ul(m))*sn(m,kdir,1) &
           +(rhor(m)*vr(m)+rhol(m)*vl(m))*sn(m,kdir,2) &
           +(rhor(m)*wr(m)+rhol(m)*wl(m))*sn(m,kdir,3)-dj1
         g2=dfxx*sn(m,kdir,1)+dfxy*sn(m,kdir,2)+dfxz*sn(m,kdir,3)-dj2
         g3=dfxy*sn(m,kdir,1)+dfyy*sn(m,kdir,2)+dfyz*sn(m,kdir,3)-dj3
         g4=dfxz*sn(m,kdir,1)+dfyz*sn(m,kdir,2)+dfzz*sn(m,kdir,3)-dj4
         g5=dfex*sn(m,kdir,1)+dfey*sn(m,kdir,2)+dfez*sn(m,kdir,3)-dj5
         u(n,1)=u(n,1)-0.5*g1
         u(n,2)=u(n,2)-0.5*g2
         u(n,3)=u(n,3)-0.5*g3
         u(n,4)=u(n,4)-0.5*g4
         u(n,5)=u(n,5)-0.5*g5
         u(n-ninc,1)=u(n-ninc,1)+0.5*g1
         u(n-ninc,2)=u(n-ninc,2)+0.5*g2
         u(n-ninc,3)=u(n-ninc,3)+0.5*g3
         u(n-ninc,4)=u(n-ninc,4)+0.5*g4
         u(n-ninc,5)=u(n-ninc,5)+0.5*g5
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
!------direction k-------------------------------------------------------
!
      if(equat(3:4).eq.'3d') then
       kdir=3
       ninc=nck
!
      do k=k1p1,k2m1
       do j=j1,j2m1
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
      do k=k1p1,k2m1
       do j=j1,j2m1
        ind1 = indc(i1,j,k)
        ind2 = indc(i2m1,j,k)
!DEC$ IVDEP
        do n=ind1,ind2
         m=n-n0c
!        vecteur normal unitaire a la face consideree
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
         el=pl(m)/(gam1*rhol(m))+0.5*q2l+pinfl/rhol(m)
         er=prr(m)/(gam1*rhor(m))+0.5*q2r+pinfl/rhor(m)
         vnl=ul(m)*nx+vl(m)*ny+wl(m)*nz
         vnr=ur(m)*nx+vr(m)*ny+wr(m)*nz
!        calcul des vitesses caracteristiques preconditionnees
         beta2l=min(max(q2l/al**2,cte*(qinf/al)**2),1.)
         beta2r=min(max(q2r/ar**2,cte*(qinf/ar)**2),1.)
         cpl=0.5*((1.+beta2l)*vnl*cnds+sqrt(((1.-beta2l)*vnl)**2+4.*beta2l*cnds**2*al**2))
         cpr=0.5*((1.+beta2r)*vnr*cnds+sqrt(((1.-beta2r)*vnr)**2+4.*beta2r*cnds**2*ar**2))
!        calcul rayon spectral
         rm=max(abs(cpl),abs(cpr))
!        evaluation du terme de dissipation (x2)
         dk1=rm*(rhor(m)-rhol(m))
         dk2=rm*(rhor(m)*ur(m)-rhol(m)*ul(m))
         dk3=rm*(rhor(m)*vr(m)-rhol(m)*vl(m))
         dk4=rm*(rhor(m)*wr(m)-rhol(m)*wl(m))
         dk5=rm*(rhor(m)*er-rhol(m)*el)
!        calcul du flux a l'interface k-1/2 (x2)
         dfxx=rhor(m)*ur(m)**2 + prr(m)-pinfl - toxx(n) &
             +rhol(m)*ul(m)**2 + pl(m)-pinfl - toxx(n-ninc)
         dfxy=rhor(m)*ur(m)*vr(m) - toxy(n) &
             +rhol(m)*ul(m)*vl(m) - toxy(n-ninc)
         dfxz=rhor(m)*ur(m)*wr(m) - toxz(n) &
             +rhol(m)*ul(m)*wl(m) - toxz(n-ninc)
         dfyy=rhor(m)*vr(m)**2 + prr(m)-pinfl - toyy(n) &
             +rhol(m)*vl(m)**2 + pl(m)-pinfl - toyy(n-ninc)
         dfyz=rhor(m)*vr(m)*wr(m) - toyz(n) &
             +rhol(m)*vl(m)*wl(m) - toyz(n-ninc)
         dfzz=rhor(m)*wr(m)**2 + prr(m)-pinfl - tozz(n) &
             +rhol(m)*wl(m)**2 + pl(m)-pinfl - tozz(n-ninc)
         dfex=(rhor(m)*er+prr(m)-pinfl)*ur(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*ul(m) &
             -(toxx(n     )*ur(m)+toxy(n     )*vr(m)+toxz(n     )*wr(m)) &
             -(toxx(n-ninc)*ul(m)+toxy(n-ninc)*vl(m)+toxz(n-ninc)*wl(m)) &
             - qcx(n) - qcx(n-ninc)
         dfey=(rhor(m)*er+prr(m)-pinfl)*vr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*vl(m) &
             -(toxy(n     )*ur(m)+toyy(n     )*vr(m)+toyz(n     )*wr(m)) &
             -(toxy(n-ninc)*ul(m)+toyy(n-ninc)*vl(m)+toyz(n-ninc)*wl(m)) &
             - qcy(n) - qcy(n-ninc)
         dfez=(rhor(m)*er+prr(m)-pinfl)*wr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*wl(m) &
             -(toxz(n     )*ur(m)+toyz(n     )*vr(m)+tozz(n     )*wr(m)) &
             -(toxz(n-ninc)*ul(m)+toyz(n-ninc)*vl(m)+tozz(n-ninc)*wl(m)) &
             - qcz(n) - qcz(n-ninc)
         h1=(rhor(m)*ur(m)+rhol(m)*ul(m))*sn(m,kdir,1) &
           +(rhor(m)*vr(m)+rhol(m)*vl(m))*sn(m,kdir,2) &
           +(rhor(m)*wr(m)+rhol(m)*wl(m))*sn(m,kdir,3)-dk1
         h2=dfxx*sn(m,kdir,1)+dfxy*sn(m,kdir,2)+dfxz*sn(m,kdir,3)-dk2
         h3=dfxy*sn(m,kdir,1)+dfyy*sn(m,kdir,2)+dfyz*sn(m,kdir,3)-dk3
         h4=dfxz*sn(m,kdir,1)+dfyz*sn(m,kdir,2)+dfzz*sn(m,kdir,3)-dk4
         h5=dfex*sn(m,kdir,1)+dfey*sn(m,kdir,2)+dfez*sn(m,kdir,3)-dk5
!
         u(n,1)=u(n,1)-0.5*h1
         u(n,2)=u(n,2)-0.5*h2
         u(n,3)=u(n,3)-0.5*h3
         u(n,4)=u(n,4)-0.5*h4
         u(n,5)=u(n,5)-0.5*h5
         u(n-ninc,1)=u(n-ninc,1)+0.5*h1
         u(n-ninc,2)=u(n-ninc,2)+0.5*h2
         u(n-ninc,3)=u(n-ninc,3)+0.5*h3
         u(n-ninc,4)=u(n-ninc,4)+0.5*h4
         u(n-ninc,5)=u(n-ninc,5)+0.5*h5
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
      if(isortie.eq.1) then
      write(6,'("===>sch_hllc: ecriture increment expli")')
      k=1
      i=80
      do j=j1,j2m1
       n=indc(i,j,k)
       m=n-n0c
       write(6,'(i4,i6,4(1pe12.4))') &
          j,n,u(n,1),u(n,2),u(n,4),u(n,5)
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
      end subroutine
end module
