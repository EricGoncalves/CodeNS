module mod_sch_weno3pond
  implicit none
contains
  subroutine sch_weno3pond( &
       lm,ityprk, &
       u,v,ff, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       equat, &
       sn,lgsnlt, &
       fxx,fyy,fzz,fxy,fxz,fyz,fex,fey,fez, &
       ps, &
       cvi,cvj,cmui1,cmui2,cmuj1,cmuj2)
!
!******************************************************************
!
!_DA  DATE_C : decembre 2004 - Eric Goncalves / LEGI
!
!     ACT
!_A   Schema WENO de Jiang&Chu, ordre 3 pondere.
!_A   Formulation 2D avec schema de Roe.
!_A   Prise en compte de la metrique du maillage.
!_A   Mapping de Henrik.
!
!*******************************************************************
!-----parameters figes----------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    implicit none
    integer          ::       i,     i1,   i1m1,   i1p1,     i2
    integer          ::    i2m1,   i2m2,     id,   iexp,   imap
    integer          ::    ind1,   ind2,isortie, ityprk,      j
    integer          ::      j1,   j1m1,   j1p1,     j2,   j2m1
    integer          ::    j2m2,     jd,      k,     k1,   k1m1
    integer          ::    k1p1,     k2,   k2m1,   k2m2,     kd
    integer          ::    kdir, lgsnlt,     lm,      m,     m1
    integer          ::       n,    n0c,     n1,    nci,    ncj
    integer          ::     nid,   nijd,   ninc,    njd
    double precision ::                   al,                  am,                am2i,                  ar,              beta11
    double precision ::               beta12,              beta21,              beta22,              beta31,              beta32
    double precision ::               beta41,              beta42,              beta51,              beta52,                 c00
    double precision ::                  c01,                 c10,                 c11,                 c20,                 c21
    double precision ::          cmui1(ip21),         cmui2(ip21),         cmuj1(ip21),         cmuj2(ip21),                cnds
    double precision ::            cvi(ip21),           cvj(ip21),                d1c1,                d1c2,                d1c3
    double precision ::                 d2c1,                d2c2,                d2c3,                d3c1,                d3c2
    double precision ::                 d3c3,                d4c1,                d4c2,                d4c3,                d5c1
    double precision ::                 d5c2,                d5c3,                 df1,                 df2,                 df3
    double precision ::                  df4,                 df5,                 dg1,                 dg2,                 dg3
    double precision ::                  dg4,                 dg5,                 eps,                  f1,                 f11
    double precision ::                  f12,                  f2,                 f21,                 f22,                  f3
    double precision ::                  f31,                 f32,                  f4,                 f41,                 f42
    double precision ::                   f5,                 f51,                 f52,                 fc1,                 fc2
    double precision ::                  fc3,                 fc4,                 fc5,           fex(ip00),           fey(ip00)
    double precision ::            fez(ip00),       ff(ip11,ip60),                 fv2,                 fv3,                 fv4
    double precision ::                  fv5,           fxx(ip00),           fxy(ip00),           fxz(ip00),           fyy(ip00)
    double precision ::            fyz(ip00),           fzz(ip00),                  g1,                 g11,                 g12
    double precision ::                  g1m,                 g1p,                  g2,                 g21,                 g22
    double precision ::                  g2m,                 g2p,                  g3,                 g31,                 g32
    double precision ::                   g4,                 g41,                 g42,                  g5,                 g51
    double precision ::                  g52,                 gc1,                 gc2,                 gc3,                 gc4
    double precision ::                  gc5,                  gd,                 gd1,                 gd2,                 gv2
    double precision ::                  gv3,                 gv4,                 gv5,                  hl,                  hm
    double precision ::                   hr,                  nx,                  ny,                  nz,                 p11
    double precision ::                  p12,                 p13,                 p14,                 p15,                 p21
    double precision ::                  p22,                 p23,                 p24,                 p25,                 p31
    double precision ::                  p32,                 p33,                 p34,                 p35,                 p41
    double precision ::                  p42,                 p43,                 p44,                 p45,                 p51
    double precision ::                  p52,                 p53,                 p54,                 p55,            ps(ip11)
    double precision ::                  q11,                 q12,                 q13,                 q14,                 q15
    double precision ::                  q1f,               q1f1m,               q1f1p,               q1f2p,                 q21
    double precision ::                  q22,                 q23,                 q24,                 q25,                 q2f
    double precision ::                q2f1m,               q2f1p,               q2f2p,                 q31,                 q32
    double precision ::                  q33,                 q34,                 q35,                 q3f,               q3f1m
    double precision ::                q3f1p,               q3f2p,                 q41,                 q42,                 q43
    double precision ::                  q44,                 q45,                 q4f,               q4f1m,               q4f1p
    double precision ::                q4f2p,                 q51,                 q52,                 q53,                 q54
    double precision ::                  q55,                 q5f,               q5f1m,               q5f1p,               q5f2p
    double precision ::            qcx(ip12),           qcy(ip12),           qcz(ip12),              rhoami,              rhoiam
    double precision ::                 rhom,               rhomi,sn(lgsnlt,nind,ndir),                  sw,                 swm
    double precision ::           toxx(ip12),          toxy(ip12),          toxz(ip12),          toyy(ip12),          toyz(ip12)
    double precision ::           tozz(ip12),        u(ip11,ip60),                  ul,                  um,                  ur
    double precision ::         v(ip11,ip60),                  v1,                  v4,                  v5,               vitm2
    double precision ::                   vl,                  vm,                  vn,                  vr,                 w11
    double precision ::                  w12,                 w13,                 w14,                 w15,                 w21
    double precision ::                  w22,                 w23,                 w24,                 w25,                  wl
    double precision ::                   wm,                  wr,                ww11,               ww11m,                ww12
    double precision ::                ww12m,                ww13,               ww13m,                ww14,               ww14m
    double precision ::                 ww15,               ww15m,                ww21,               ww21m,                ww22
    double precision ::                ww22m,                ww23,               ww23m,                ww24,               ww24m
    double precision ::                 ww25,               ww25m
!
!-------------------------------------------------------------------
!
    character(len=7 ) :: equat




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
    i2m2=i2-2
    j2m2=j2-2
    k2m2=k2-2
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
!
    nci = inc(1,0,0)
    ncj = inc(0,1,0)
!
!     type de maillage=maillage (x,y)
!     activation des sorties
    isortie=0
!     activation mapping
    imap=1
!     epsilon petit
    eps=1.e-40
!      eps=1.e-6
!
!-----calcul des densites de flux convectifs ----------------------
!
    ind1 = indc(i1m1,j1m1,k1  )
    ind2 = indc(i2  ,j2  ,k2m1)
!$OMP SIMD
    do n=ind1,ind2
       m=n-n0c
       u(n,1)=0.
       u(n,2)=0.
       u(n,3)=0.
       u(n,4)=0.
       u(n,5)=0.
       fxx(m)=v(n,2)*(v(n,2)/v(n,1))+ps(n)-pinfl
       fxy(m)=v(n,3)*(v(n,2)/v(n,1))
       fxz(m)=v(n,4)*(v(n,2)/v(n,1))
       fyy(m)=v(n,3)*(v(n,3)/v(n,1))+ps(n)-pinfl
       fyz(m)=v(n,4)*(v(n,3)/v(n,1))
       fzz(m)=v(n,4)*(v(n,4)/v(n,1))+ps(n)-pinfl
       fex(m)=(v(n,5)+ps(n)-pinfl)*v(n,2)/v(n,1)
       fey(m)=(v(n,5)+ps(n)-pinfl)*v(n,3)/v(n,1)
       fez(m)=(v(n,5)+ps(n)-pinfl)*v(n,4)/v(n,1)
    enddo
!
!*****************************************************************
!
!  Calcul du flux numerique par direction suivant les etapessuccessives :
!    1) evaluation des matrices de passage et des valeurs
!       propres associees a la matrice jacobienne A
!       (ces quantites sont evaluees a l'etat moyen de Roe)
!    2) calculs des flux associes aux variables caracteristiques
!    3) application de la procedure de reconstruction ENO scalaire
!       a chaque composante
!    4) calculs des flux convectifs a partir des flux reconstruits
!    5) ajout des flux visqueux evalues avec schema centre
!
!********************************************************************

    if(imap.eq.0) then
!
!-----direction i-----------------------------------------
!
       kdir=1
       ninc=nci
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m2,j,k)
             do n=ind1,ind2
                m=n-n0c
                m1=m+ninc
                n1=n+ninc
!        vecteur normal unitaire a la face consideree (face i+1/2)
                cnds=sqrt(sn(m1,kdir,1)*sn(m1,kdir,1)+ &
                     sn(m1,kdir,2)*sn(m1,kdir,2)+ &
                     sn(m1,kdir,3)*sn(m1,kdir,3))
                nx=sn(m1,kdir,1)/cnds
                ny=sn(m1,kdir,2)/cnds
                nz=sn(m1,kdir,3)/cnds
!        calcul des etats gauche et droit
                ul=v(n,2)/v(n,1)
                vl=v(n,3)/v(n,1)
                wl=v(n,4)/v(n,1)
                ur=v(n1,2)/v(n1,1)
                vr=v(n1,3)/v(n1,1)
                wr=v(n1,4)/v(n1,1)
                al=sqrt(gam*ps(n )/v(n,1))
                ar=sqrt(gam*ps(n1)/v(n1,1))
                hl=al*al/gam1+0.5*(ul**2+vl**2+wl**2)
                hr=ar*ar/gam1+0.5*(ur**2+vr**2+wr**2)
!        calcul de etat moyen de Roe
                gd=sqrt(v(n1,1)/v(n,1))
                gd1=1./(1.+gd)
                gd2=gd*gd1
                rhom=sqrt(v(n,1)*v(n1,1))
                rhomi=1./rhom
                um=gd1*ul+gd2*ur
                vm=gd1*vl+gd2*vr
                wm=gd1*wl+gd2*wr
                hm=gd1*hl+gd2*hr
                vitm2=0.5*(um**2+vm**2+wm**2)
                am=sqrt(abs(gam1*(hm-vitm2)))
                am2i=1./(am*am)
                vn=um*nx+vm*ny+wm*nz
                rhoiam=rhom/am
                rhoami=am2i/rhoiam
!        valeurs propres de la matrice jacobienne des flux
                v1=vn
                v4=vn+am
                v5=vn-am
!        calcul des matrices de passage a gauche Q et a droite P
                q11=(1.-gam1*vitm2*am2i)*nx-(vm*nz-wm*ny)*rhomi
                q12=gam1*um*nx*am2i
                q13=gam1*vm*nx*am2i+nz*rhomi
                q14=gam1*wm*nx*am2i-ny*rhomi
                q15=-gam1*nx*am2i
                q21=(1.-gam1*vitm2*am2i)*ny-(wm*nx-um*nz)*rhomi
                q22=gam1*um*ny*am2i-nz*rhomi
                q23=gam1*vm*ny*am2i
                q24=gam1*wm*ny*am2i+nx*rhomi
                q25=-gam1*ny*am2i
                q31=(1.-gam1*vitm2*am2i)*nz-(um*ny-vm*nx)*rhomi
                q32=gam1*um*nz*am2i+ny*rhomi
                q33=gam1*vm*nz*am2i-nx*rhomi
                q34=gam1*wm*nz*am2i
                q35=-gam1*nz*am2i
                q41=gam1*vitm2*rhoami-vn*rhomi
                q42=nx*rhomi-gam1*um*rhoami
                q43=ny*rhomi-gam1*vm*rhoami
                q44=nz*rhomi-gam1*wm*rhoami
                q45=gam1*rhoami
                q51=gam1*vitm2*rhoami+vn*rhomi
                q52=-nx*rhomi-gam1*um*rhoami
                q53=-ny*rhomi-gam1*vm*rhoami
                q54=-nz*rhomi-gam1*wm*rhoami
                q55=gam1*rhoami
!
                p11=nx
                p12=ny
                p13=nz
                p14=0.5*rhoiam
                p15=0.5*rhoiam
                p21=um*nx
                p22=um*ny-rhom*nz
                p23=um*nz+rhom*ny
                p24=0.5*rhoiam*(um+nx*am)
                p25=0.5*rhoiam*(um-nx*am)
                p31=vm*nx+rhom*nz
                p32=vm*ny
                p33=vm*nz-rhom*nx
                p34=0.5*rhoiam*(vm+ny*am)
                p35=0.5*rhoiam*(vm-ny*am)
                p41=wm*nx-rhom*ny
                p42=wm*ny+rhom*nx
                p43=wm*nz
                p44=0.5*rhoiam*(wm+nz*am)
                p45=0.5*rhoiam*(wm-nz*am)
                p51=vitm2*nx+rhom*(vm*nz-wm*ny)
                p52=vitm2*ny+rhom*(wm*nx-um*nz)
                p53=vitm2*nz+rhom*(um*ny-vm*nx)
                p54=0.5*rhoiam*(hm+am*vn)
                p55=0.5*rhoiam*(hm-am*vn)
!        produit de Q avec les flux Euler aux points (i-1) a (i+2)
                q1f1m=q11*v(n-ninc  ,2)+q12*fxx(m-ninc)  +q13*fxy(m-ninc) &
                     +q14*fxz(m-ninc)  +q15*fex(m-ninc)
                q1f  =q11*v(n       ,2)+q12*fxx(m)       +q13*fxy(m) &
                     +q14*fxz(m)       +q15*fex(m)
                q1f1p=q11*v(n+ninc  ,2)+q12*fxx(m+ninc)  +q13*fxy(m+ninc) &
                     +q14*fxz(m+ninc)  +q15*fex(m+ninc)
                q1f2p=q11*v(n+2*ninc,2)+q12*fxx(m+2*ninc)+q13*fxy(m+2*ninc) &
                     +q14*fxz(m+2*ninc)+q15*fex(m+2*ninc)
!
                q2f1m=q21*v(n-ninc  ,2)+q22*fxx(m-ninc)  +q23*fxy(m-ninc) &
                     +q24*fxz(m-ninc)  +q25*fex(m-ninc)
                q2f  =q21*v(n       ,2)+q22*fxx(m)       +q23*fxy(m) &
                     +q24*fxz(m)       +q25*fex(m)
                q2f1p=q21*v(n+ninc  ,2)+q22*fxx(m+ninc)  +q23*fxy(m+ninc) &
                     +q24*fxz(m+ninc)  +q25*fex(m+ninc)
                q2f2p=q21*v(n+2*ninc,2)+q22*fxx(m+2*ninc)+q23*fxy(m+2*ninc) &
                     +q24*fxz(m+2*ninc)+q25*fex(m+2*ninc)
!
                q3f1m=q31*v(n-ninc  ,2)+q32*fxx(m-ninc)  +q33*fxy(m-ninc) &
                     +q34*fxz(m-ninc)  +q35*fex(m-ninc)
                q3f  =q31*v(n       ,2)+q32*fxx(m)       +q33*fxy(m) &
                     +q34*fxz(m)       +q35*fex(m)
                q3f1p=q31*v(n+ninc  ,2)+q32*fxx(m+ninc)  +q33*fxy(m+ninc) &
                     +q34*fxz(m+ninc)  +q35*fex(m+ninc)
                q3f2p=q31*v(n+2*ninc,2)+q32*fxx(m+2*ninc)+q33*fxy(m+2*ninc) &
                     +q34*fxz(m+2*ninc)+q35*fex(m+2*ninc)
!
                q4f1m=q41*v(n-ninc  ,2)+q42*fxx(m-ninc)  +q43*fxy(m-ninc) &
                     +q44*fxz(m-ninc)  +q45*fex(m-ninc)
                q4f  =q41*v(n       ,2)+q42*fxx(m)       +q43*fxy(m) &
                     +q44*fxz(m)       +q45*fex(m)
                q4f1p=q41*v(n+ninc  ,2)+q42*fxx(m+ninc)  +q43*fxy(m+ninc) &
                     +q44*fxz(m+ninc)  +q45*fex(m+ninc)
                q4f2p=q41*v(n+2*ninc,2)+q42*fxx(m+2*ninc)+q43*fxy(m+2*ninc) &
                     +q44*fxz(m+2*ninc)+q45*fex(m+2*ninc)
!
                q5f1m=q51*v(n-ninc  ,2)+q52*fxx(m-ninc)  +q53*fxy(m-ninc) &
                     +q54*fxz(m-ninc)  +q55*fex(m-ninc)
                q5f  =q51*v(n       ,2)+q52*fxx(m)       +q53*fxy(m) &
                     +q54*fxz(m)       +q55*fex(m)
                q5f1p=q51*v(n+ninc  ,2)+q52*fxx(m+ninc)  +q53*fxy(m+ninc) &
                     +q54*fxz(m+ninc)  +q55*fex(m+ninc)
                q5f2p=q51*v(n+2*ninc,2)+q52*fxx(m+2*ninc)+q53*fxy(m+2*ninc) &
                     +q54*fxz(m+2*ninc)+q55*fex(m+2*ninc)
!        coefficients Crj en maillage irregulier
                c00=cmui2(m1)*0.5
                c01=cmui1(m1)*0.5
                c10=-cmui2(m)*0.5
                c11=1.-c10
                c21=-cmui1(m+2*ninc)*0.5
                c20=1.-c21
!        calcul des flux d'ordre 2 sur les 2 stencils
                f11=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q1f  *c00 +q1f1p*c01)
                f12=0.5*(1.+sign(1.,v1))*(q1f  *c00 +q1f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q1f1p*c11 +q1f2p*c10)
!
                f21=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q2f  *c00 +q2f1p*c01)
                f22=0.5*(1.+sign(1.,v1))*(q2f  *c00 +q2f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q2f1p*c11 +q2f2p*c10)
!
                f31=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q3f  *c00 +q3f1p*c01)
                f32=0.5*(1.+sign(1.,v1))*(q3f  *c00 +q3f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q3f1p*c11 +q3f2p*c10)
!
                f41=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11) &
                     +0.5*(1.-sign(1.,v4))*(q4f  *c00 +q4f1p*c01)
                f42=0.5*(1.+sign(1.,v4))*(q4f  *c00 +q4f1p*c01) &
                     +0.5*(1.-sign(1.,v4))*(q4f1p*c11 +q4f2p*c10)
!
                f51=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11) &
                     +0.5*(1.-sign(1.,v5))*(q5f  *c00 +q5f1p*c01)
                f52=0.5*(1.+sign(1.,v5))*(q5f  *c00 +q5f1p*c01) &
                     +0.5*(1.-sign(1.,v5))*(q5f1p*c11 +q5f2p*c10)
!        calcul des senseurs beta (au carre)
!         d1c1=(q1f-q1f1m)**2    !maillage regulier
!         d1c2=(q1f1p-q1f)**2
!         d1c3=(q1f2p-q1f1p)**2
                d1c1=4.*((q1f-q1f1m)*c10)**2
                d1c2=4.*((q1f1p-q1f)*c01)**2
                d1c3=4.*((q1f2p-q1f1p)*c01*c21/c00)**2
                d2c1=4.*((q2f-q2f1m)*c10)**2
                d2c2=4.*((q2f1p-q2f)*c01)**2
                d2c3=4.*((q2f2p-q2f1p)*c01*c21/c00)**2
                d3c1=4.*((q3f-q3f1m)*c10)**2
                d3c2=4.*((q3f1p-q3f)*c01)**2
                d3c3=4.*((q3f2p-q3f1p)*c01*c21/c00)**2
                d4c1=4.*((q4f-q4f1m)*c10)**2
                d4c2=4.*((q4f1p-q4f)*c01)**2
                d4c3=4.*((q4f2p-q4f1p)*c01*c21/c00)**2
                d5c1=4.*((q5f-q5f1m)*c10)**2
                d5c2=4.*((q5f1p-q5f)*c01)**2
                d5c3=4.*((q5f2p-q5f1p)*c01*c21/c00)**2
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*d1c1 &
                     +0.5*(1.-sign(1.,v1))*d1c2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*d1c2 &
                     +0.5*(1.-sign(1.,v1))*d1c3+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*d2c1 &
                     +0.5*(1.-sign(1.,v1))*d2c2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*d2c2 &
                     +0.5*(1.-sign(1.,v1))*d2c3+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*d3c1 &
                     +0.5*(1.-sign(1.,v1))*d3c2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*d3c2 &
                     +0.5*(1.-sign(1.,v1))*d3c3+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*d4c1 &
                     +0.5*(1.-sign(1.,v4))*d4c2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*d4c2 &
                     +0.5*(1.-sign(1.,v4))*d4c3+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*d5c1 &
                     +0.5*(1.-sign(1.,v5))*d5c2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*d5c2 &
                     +0.5*(1.-sign(1.,v5))*d5c3+eps)**iexp
!        coefficients gamma en maillage irregulier
                g1p=cmui2(m1)*cvi(m1)/(cmui1(m)*cvi(m)+cmui2(m)*cvi(m)+cmui2(m1)*cvi(m1))
                g2p=1.-g1p
                g2m=cmui1(m1)*cvi(m1)/(cmui1(m1)*cvi(m1)+cmui2(m1)*cvi(m1)+cmui2(m+2*ninc)*cvi(m+2*ninc))
                g1m=1.-g2m
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
!
                ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
!
                ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
!
                ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                     +0.5*(1.-sign(1.,v4))*(g1m/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                     +0.5*(1.-sign(1.,v4))*(g2m/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
!
                ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                     +0.5*(1.-sign(1.,v5))*(g1m/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                     +0.5*(1.-sign(1.,v5))*(g2m/beta52)
                sw=ww15+ww25
                w15=ww15/sw
                w25=ww25/sw
!        calcul des flux convectifs projetes
                fc1=w11*f11+w21*f12
                fc2=w12*f21+w22*f22
                fc3=w13*f31+w23*f32
                fc4=w14*f41+w24*f42
                fc5=w15*f51+w25*f52
!        produit avec matrice P pour retour dans l'espace physique
                f1=fc1*p11+fc2*p12+fc3*p13+fc4*p14+fc5*p15
                f2=fc1*p21+fc2*p22+fc3*p23+fc4*p24+fc5*p25
                f3=fc1*p31+fc2*p32+fc3*p33+fc4*p34+fc5*p35
                f5=fc1*p51+fc2*p52+fc3*p53+fc4*p54+fc5*p55
!-------------------------------------------------------------------
!        produit de Q avec les flux Euler aux points (i-1) a (i+2)
                q1f1m=q11*v(n-ninc  ,3)+q12*fxy(m-ninc)  +q13*fyy(m-ninc) &
                     +q14*fyz(m-ninc)  +q15*fey(m-ninc)
                q1f  =q11*v(n       ,3)+q12*fxy(m)       +q13*fyy(m) &
                     +q14*fyz(m)       +q15*fey(m)
                q1f1p=q11*v(n+ninc  ,3)+q12*fxy(m+ninc)  +q13*fyy(m+ninc) &
                     +q14*fyz(m+ninc)  +q15*fey(m+ninc)
                q1f2p=q11*v(n+2*ninc,3)+q12*fxy(m+2*ninc)+q13*fyy(m+2*ninc) &
                     +q14*fyz(m+2*ninc)+q15*fey(m+2*ninc)
!
                q2f1m=q21*v(n-ninc  ,3)+q22*fxy(m-ninc)  +q23*fyy(m-ninc) &
                     +q24*fyz(m-ninc)  +q25*fey(m-ninc)
                q2f  =q21*v(n       ,3)+q22*fxy(m)       +q23*fyy(m) &
                     +q24*fyz(m)       +q25*fey(m)
                q2f1p=q21*v(n+ninc  ,3)+q22*fxy(m+ninc)  +q23*fyy(m+ninc) &
                     +q24*fyz(m+ninc)  +q25*fey(m+ninc)
                q2f2p=q21*v(n+2*ninc,3)+q22*fxy(m+2*ninc)+q23*fyy(m+2*ninc) &
                     +q24*fyz(m+2*ninc)+q25*fey(m+2*ninc)
!
                q3f1m=q31*v(n-ninc  ,3)+q32*fxy(m-ninc)  +q33*fyy(m-ninc) &
                     +q34*fyz(m-ninc)  +q35*fey(m-ninc)
                q3f  =q31*v(n       ,3)+q32*fxy(m)       +q33*fyy(m) &
                     +q34*fyz(m)       +q35*fey(m)
                q3f1p=q31*v(n+ninc  ,3)+q32*fxy(m+ninc)  +q33*fyy(m+ninc) &
                     +q34*fyz(m+ninc)  +q35*fey(m+ninc)
                q3f2p=q31*v(n+2*ninc,3)+q32*fxy(m+2*ninc)+q33*fyy(m+2*ninc) &
                     +q34*fyz(m+2*ninc)+q35*fey(m+2*ninc)
!
                q4f1m=q41*v(n-ninc  ,3)+q42*fxy(m-ninc)  +q43*fyy(m-ninc) &
                     +q44*fyz(m-ninc)  +q45*fey(m-ninc)
                q4f  =q41*v(n       ,3)+q42*fxy(m)       +q43*fyy(m) &
                     +q44*fyz(m)       +q45*fey(m)
                q4f1p=q41*v(n+ninc  ,3)+q42*fxy(m+ninc)  +q43*fyy(m+ninc) &
                     +q44*fyz(m+ninc)  +q45*fey(m+ninc)
                q4f2p=q41*v(n+2*ninc,3)+q42*fxy(m+2*ninc)+q43*fyy(m+2*ninc) &
                     +q44*fyz(m+2*ninc)+q45*fey(m+2*ninc)
!
                q5f1m=q51*v(n-ninc  ,3)+q52*fxy(m-ninc)  +q53*fyy(m-ninc) &
                     +q54*fyz(m-ninc)  +q55*fey(m-ninc)
                q5f  =q51*v(n       ,3)+q52*fxy(m)       +q53*fyy(m) &
                     +q54*fyz(m)       +q55*fey(m)
                q5f1p=q51*v(n+ninc  ,3)+q52*fxy(m+ninc)  +q53*fyy(m+ninc) &
                     +q54*fyz(m+ninc)  +q55*fey(m+ninc)
                q5f2p=q51*v(n+2*ninc,3)+q52*fxy(m+2*ninc)+q53*fyy(m+2*ninc) &
                     +q54*fyz(m+2*ninc)+q55*fey(m+2*ninc)
!        calcul des flux d'ordre 2 sur les 2 stencils
                g11=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q1f  *c00 +q1f1p*c01)
                g12=0.5*(1.+sign(1.,v1))*(q1f  *c00 +q1f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q1f1p*c11 +q1f2p*c10)
!
                g21=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q2f  *c00 +q2f1p*c01)
                g22=0.5*(1.+sign(1.,v1))*(q2f  *c00 +q2f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q2f1p*c11 +q2f2p*c10)
!
                g31=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q3f  *c00 +q3f1p*c01)
                g32=0.5*(1.+sign(1.,v1))*(q3f  *c00 +q3f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q3f1p*c11 +q3f2p*c10)
!
                g41=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11) &
                     +0.5*(1.-sign(1.,v4))*(q4f  *c00 +q4f1p*c01)
                g42=0.5*(1.+sign(1.,v4))*(q4f  *c00 +q4f1p*c01) &
                     +0.5*(1.-sign(1.,v4))*(q4f1p*c11 +q4f2p*c10)
!
                g51=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11) &
                     +0.5*(1.-sign(1.,v5))*(q5f  *c00 +q5f1p*c01)
                g52=0.5*(1.+sign(1.,v5))*(q5f  *c00 +q5f1p*c01) &
                     +0.5*(1.-sign(1.,v5))*(q5f1p*c11 +q5f2p*c10)
!        calcul des senseurs beta (au carre)
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*d1c1 &
                     +0.5*(1.-sign(1.,v1))*d1c2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*d1c2 &
                     +0.5*(1.-sign(1.,v1))*d1c3+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*d2c1 &
                     +0.5*(1.-sign(1.,v1))*d2c2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*d2c2 &
                     +0.5*(1.-sign(1.,v1))*d2c3+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*d3c1 &
                     +0.5*(1.-sign(1.,v1))*d3c2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*d3c2 &
                     +0.5*(1.-sign(1.,v1))*d3c3+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*d4c1 &
                     +0.5*(1.-sign(1.,v4))*d4c2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*d4c2 &
                     +0.5*(1.-sign(1.,v4))*d4c3+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*d5c1 &
                     +0.5*(1.-sign(1.,v5))*d5c2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*d5c2 &
                     +0.5*(1.-sign(1.,v5))*d5c3+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
!
                ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
!
                ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
!
                ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                     +0.5*(1.-sign(1.,v4))*(g1m/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                     +0.5*(1.-sign(1.,v4))*(g2m/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
!
                ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                     +0.5*(1.-sign(1.,v5))*(g1m/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                     +0.5*(1.-sign(1.,v5))*(g2m/beta52)
                sw=ww15+ww25
                w15=ww15/sw
                w25=ww25/sw
!        calcul des flux convectifs projetes
                gc1=w11*g11+w21*g12
                gc2=w12*g21+w22*g22
                gc3=w13*g31+w23*g32
                gc4=w14*g41+w24*g42
                gc5=w15*g51+w25*g52
!        produit avec matrice P pour retour dans l'espace physique
                g1=gc1*p11+gc2*p12+gc3*p13+gc4*p14+gc5*p15
                g2=gc1*p21+gc2*p22+gc3*p23+gc4*p24+gc5*p25
                g3=gc1*p31+gc2*p32+gc3*p33+gc4*p34+gc5*p35
                g5=gc1*p51+gc2*p52+gc3*p53+gc4*p54+gc5*p55
!        calcul du flux numerique et bilan de flux
                df1=f1*sn(m1,kdir,1)+g1*sn(m1,kdir,2)
                df2=f2*sn(m1,kdir,1)+g2*sn(m1,kdir,2)
                df3=f3*sn(m1,kdir,1)+g3*sn(m1,kdir,2)
                df5=f5*sn(m1,kdir,1)+g5*sn(m1,kdir,2)
!        calcul des flux visqueux (multiplies par -2)
                fv2=(cmui2(m1)*toxx(n)+cmui1(m1)*toxx(n1))*sn(m1,kdir,1) &
                     +(cmui2(m1)*toxy(n)+cmui1(m1)*toxy(n1))*sn(m1,kdir,2)
                fv3=(cmui2(m1)*toxy(n)+cmui1(m1)*toxy(n1))*sn(m1,kdir,1) &
                     +(cmui2(m1)*toyy(n)+cmui1(m1)*toyy(n1))*sn(m1,kdir,2)
                fv5=(cmui2(m1)*(toxx(n )*ul+toxy(n )*vl+qcx(n )) &
                     +cmui1(m1)*(toxx(n1)*ur+toxy(n1)*vr+qcx(n1)))*sn(m1,kdir,1) &
                     +(cmui2(m1)*(toxy(n )*ul+toyy(n )*vl+qcy(n )) &
                     +cmui1(m1)*(toxy(n1)*ur+toyy(n1)*vr+qcy(n1)))*sn(m1,kdir,2)
                u(n1,1)=u(n1,1)-df1
                u(n1,2)=u(n1,2)-df2+0.5*fv2
                u(n1,3)=u(n1,3)-df3+0.5*fv3
                u(n1,5)=u(n1,5)-df5+0.5*fv5
                u(n,1)=u(n,1)+df1
                u(n,2)=u(n,2)+df2-0.5*fv2
                u(n,3)=u(n,3)+df3-0.5*fv3
                u(n,5)=u(n,5)+df5-0.5*fv5
             enddo
          enddo
       enddo
!
!------direction j----------------------------------------------
!
       kdir=2
       ninc=ncj
!
       do k=k1,k2m1
          do j=j1,j2m2
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                m1=m+ninc
                n1=n+ninc
!        vecteur normal unitaire a la face consideree (face j+1/2)
                cnds=sqrt(sn(m1,kdir,1)*sn(m1,kdir,1)+ &
                     sn(m1,kdir,2)*sn(m1,kdir,2)+ &
                     sn(m1,kdir,3)*sn(m1,kdir,3))
                nx=sn(m1,kdir,1)/cnds
                ny=sn(m1,kdir,2)/cnds
                nz=sn(m1,kdir,3)/cnds
!        calcul des etats gauche et droit
                ul=v(n,2)/v(n,1)
                vl=v(n,3)/v(n,1)
                wl=v(n,4)/v(n,1)
                ur=v(n1,2)/v(n1,1)
                vr=v(n1,3)/v(n1,1)
                wr=v(n1,4)/v(n1,1)
                al=sqrt(gam*ps(n )/v(n,1))
                ar=sqrt(gam*ps(n1)/v(n1,1))
                hl=al*al/gam1+0.5*(ul**2+vl**2+wl**2)
                hr=ar*ar/gam1+0.5*(ur**2+vr**2+wr**2)
!        calcul des etats moyens de Roe
                gd=sqrt(v(n1,1)/v(n,1))
                gd1=1./(1.+gd)
                gd2=gd*gd1
                rhom=sqrt(v(n,1)*v(n1,1))
                rhomi=1./rhom
                um=gd1*ul+gd2*ur
                vm=gd1*vl+gd2*vr
                wm=gd1*wl+gd2*wr
                hm=gd1*hl+gd2*hr
                vitm2=0.5*(um**2+vm**2+wm**2)
                am=sqrt(abs(gam1*(hm-vitm2)))
                am2i=1./(am*am)
                vn=um*nx+vm*ny+wm*nz
                rhoiam=rhom/am
                rhoami=am2i/rhoiam
!        valeurs propres
                v1=vn
                v4=vn+am
                v5=vn-am
!        calcul des matrices de passage a gauche Q et a droite P
                q11=(1.-gam1*vitm2*am2i)*nx-(vm*nz-wm*ny)*rhomi
                q12=gam1*um*nx*am2i
                q13=gam1*vm*nx*am2i+nz*rhomi
                q14=gam1*wm*nx*am2i-ny*rhomi
                q15=-gam1*nx*am2i
                q21=(1.-gam1*vitm2*am2i)*ny-(wm*nx-um*nz)*rhomi
                q22=gam1*um*ny*am2i-nz*rhomi
                q23=gam1*vm*ny*am2i
                q24=gam1*wm*ny*am2i+nx*rhomi
                q25=-gam1*ny*am2i
                q31=(1.-gam1*vitm2*am2i)*nz-(um*ny-vm*nx)*rhomi
                q32=gam1*um*nz*am2i+ny*rhomi
                q33=gam1*vm*nz*am2i-nx*rhomi
                q34=gam1*wm*nz*am2i
                q35=-gam1*nz*am2i
                q41=gam1*vitm2*rhoami-vn*rhomi
                q42=nx*rhomi-gam1*um*rhoami
                q43=ny*rhomi-gam1*vm*rhoami
                q44=nz*rhomi-gam1*wm*rhoami
                q45=gam1*rhoami
                q51=gam1*vitm2*rhoami+vn*rhomi
                q52=-nx*rhomi-gam1*um*rhoami
                q53=-ny*rhomi-gam1*vm*rhoami
                q54=-nz*rhomi-gam1*wm*rhoami
                q55=gam1*rhoami
!
                p11=nx
                p12=ny
                p13=nz
                p14=0.5*rhoiam
                p15=0.5*rhoiam
                p21=um*nx
                p22=um*ny-rhom*nz
                p23=um*nz+rhom*ny
                p24=0.5*rhoiam*(um+nx*am)
                p25=0.5*rhoiam*(um-nx*am)
                p31=vm*nx+rhom*nz
                p32=vm*ny
                p33=vm*nz-rhom*nx
                p34=0.5*rhoiam*(vm+ny*am)
                p35=0.5*rhoiam*(vm-ny*am)
                p41=wm*nx-rhom*ny
                p42=wm*ny+rhom*nx
                p43=wm*nz
                p44=0.5*rhoiam*(wm+nz*am)
                p45=0.5*rhoiam*(wm-nz*am)
                p51=vitm2*nx+rhom*(vm*nz-wm*ny)
                p52=vitm2*ny+rhom*(wm*nx-um*nz)
                p53=vitm2*nz+rhom*(um*ny-vm*nx)
                p54=0.5*rhoiam*(hm+am*vn)
                p55=0.5*rhoiam*(hm-am*vn)
!        produit de Q avec les flux Euler aux points (i-1) a (i+2)
                q1f1m=q11*v(n-ninc  ,2)+q12*fxx(m-ninc)  +q13*fxy(m-ninc) &
                     +q14*fxz(m-ninc)  +q15*fex(m-ninc)
                q1f  =q11*v(n       ,2)+q12*fxx(m)       +q13*fxy(m) &
                     +q14*fxz(m)       +q15*fex(m)
                q1f1p=q11*v(n+ninc  ,2)+q12*fxx(m+ninc)  +q13*fxy(m+ninc) &
                     +q14*fxz(m+ninc)  +q15*fex(m+ninc)
                q1f2p=q11*v(n+2*ninc,2)+q12*fxx(m+2*ninc)+q13*fxy(m+2*ninc) &
                     +q14*fxz(m+2*ninc)+q15*fex(m+2*ninc)
!
                q2f1m=q21*v(n-ninc  ,2)+q22*fxx(m-ninc)  +q23*fxy(m-ninc) &
                     +q24*fxz(m-ninc)  +q25*fex(m-ninc)
                q2f  =q21*v(n       ,2)+q22*fxx(m)       +q23*fxy(m) &
                     +q24*fxz(m)       +q25*fex(m)
                q2f1p=q21*v(n+ninc  ,2)+q22*fxx(m+ninc)  +q23*fxy(m+ninc) &
                     +q24*fxz(m+ninc)  +q25*fex(m+ninc)
                q2f2p=q21*v(n+2*ninc,2)+q22*fxx(m+2*ninc)+q23*fxy(m+2*ninc) &
                     +q24*fxz(m+2*ninc)+q25*fex(m+2*ninc)
!
                q3f1m=q31*v(n-ninc  ,2)+q32*fxx(m-ninc)  +q33*fxy(m-ninc) &
                     +q34*fxz(m-ninc)  +q35*fex(m-ninc)
                q3f  =q31*v(n       ,2)+q32*fxx(m)       +q33*fxy(m) &
                     +q34*fxz(m)       +q35*fex(m)
                q3f1p=q31*v(n+ninc  ,2)+q32*fxx(m+ninc)  +q33*fxy(m+ninc) &
                     +q34*fxz(m+ninc)  +q35*fex(m+ninc)
                q3f2p=q31*v(n+2*ninc,2)+q32*fxx(m+2*ninc)+q33*fxy(m+2*ninc) &
                     +q34*fxz(m+2*ninc)+q35*fex(m+2*ninc)
!
                q4f1m=q41*v(n-ninc  ,2)+q42*fxx(m-ninc)  +q43*fxy(m-ninc) &
                     +q44*fxz(m-ninc)  +q45*fex(m-ninc)
                q4f  =q41*v(n       ,2)+q42*fxx(m)       +q43*fxy(m) &
                     +q44*fxz(m)       +q45*fex(m)
                q4f1p=q41*v(n+ninc  ,2)+q42*fxx(m+ninc)  +q43*fxy(m+ninc) &
                     +q44*fxz(m+ninc)  +q45*fex(m+ninc)
                q4f2p=q41*v(n+2*ninc,2)+q42*fxx(m+2*ninc)+q43*fxy(m+2*ninc) &
                     +q44*fxz(m+2*ninc)+q45*fex(m+2*ninc)
!
                q5f1m=q51*v(n-ninc  ,2)+q52*fxx(m-ninc)  +q53*fxy(m-ninc) &
                     +q54*fxz(m-ninc)  +q55*fex(m-ninc)
                q5f  =q51*v(n       ,2)+q52*fxx(m)       +q53*fxy(m) &
                     +q54*fxz(m)       +q55*fex(m)
                q5f1p=q51*v(n+ninc  ,2)+q52*fxx(m+ninc)  +q53*fxy(m+ninc) &
                     +q54*fxz(m+ninc)  +q55*fex(m+ninc)
                q5f2p=q51*v(n+2*ninc,2)+q52*fxx(m+2*ninc)+q53*fxy(m+2*ninc) &
                     +q54*fxz(m+2*ninc)+q55*fex(m+2*ninc)
!        coefficients Crj en maillage irregulier
                c00=cmuj2(m1)*0.5
                c01=cmuj1(m1)*0.5
                c10=-cmuj2(m)*0.5
                c11=1.-c10
                c21=-cmuj1(m+2*ninc)*0.5
                c20=1.-c21
!        calcul des flux d'ordre 2 sur les 2 stencils
                f11=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q1f  *c00 +q1f1p*c01)
                f12=0.5*(1.+sign(1.,v1))*(q1f  *c00 +q1f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q1f1p*c11 +q1f2p*c10)
!
                f21=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q2f  *c00 +q2f1p*c01)
                f22=0.5*(1.+sign(1.,v1))*(q2f  *c00 +q2f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q2f1p*c11 +q2f2p*c10)
!
                f31=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q3f  *c00 +q3f1p*c01)
                f32=0.5*(1.+sign(1.,v1))*(q3f  *c00 +q3f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q3f1p*c11 +q3f2p*c10)
!
                f41=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11) &
                     +0.5*(1.-sign(1.,v4))*(q4f  *c00 +q4f1p*c01)
                f42=0.5*(1.+sign(1.,v4))*(q4f  *c00 +q4f1p*c01) &
                     +0.5*(1.-sign(1.,v4))*(q4f1p*c11 +q4f2p*c10)
!
                f51=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11) &
                     +0.5*(1.-sign(1.,v5))*(q5f  *c00 +q5f1p*c01)
                f52=0.5*(1.+sign(1.,v5))*(q5f  *c00 +q5f1p*c01) &
                     +0.5*(1.-sign(1.,v5))*(q5f1p*c11 +q5f2p*c10)
!        calcul des senseurs beta (au carre)
!         d1c1=(q1f-q1f1m)**2    !maillage regulier
!         d1c2=(q1f1p-q1f)**2
!         d1c3=(q1f2p-q1f1p)**2
                d1c1=4.*((q1f-q1f1m)*c10)**2
                d1c2=4.*((q1f1p-q1f)*c01)**2
                d1c3=4.*((q1f2p-q1f1p)*c01*c21/c00)**2
                d2c1=4.*((q2f-q2f1m)*c10)**2
                d2c2=4.*((q2f1p-q2f)*c01)**2
                d2c3=4.*((q2f2p-q2f1p)*c01*c21/c00)**2
                d3c1=4.*((q3f-q3f1m)*c10)**2
                d3c2=4.*((q3f1p-q3f)*c01)**2
                d3c3=4.*((q3f2p-q3f1p)*c01*c21/c00)**2
                d4c1=4.*((q4f-q4f1m)*c10)**2
                d4c2=4.*((q4f1p-q4f)*c01)**2
                d4c3=4.*((q4f2p-q4f1p)*c01*c21/c00)**2
                d5c1=4.*((q5f-q5f1m)*c10)**2
                d5c2=4.*((q5f1p-q5f)*c01)**2
                d5c3=4.*((q5f2p-q5f1p)*c01*c21/c00)**2
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*d1c1 &
                     +0.5*(1.-sign(1.,v1))*d1c2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*d1c2 &
                     +0.5*(1.-sign(1.,v1))*d1c3+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*d2c1 &
                     +0.5*(1.-sign(1.,v1))*d2c2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*d2c2 &
                     +0.5*(1.-sign(1.,v1))*d2c3+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*d3c1 &
                     +0.5*(1.-sign(1.,v1))*d3c2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*d3c2 &
                     +0.5*(1.-sign(1.,v1))*d3c3+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*d4c1 &
                     +0.5*(1.-sign(1.,v4))*d4c2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*d4c2 &
                     +0.5*(1.-sign(1.,v4))*d4c3+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*d5c1 &
                     +0.5*(1.-sign(1.,v5))*d5c2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*d5c2 &
                     +0.5*(1.-sign(1.,v5))*d5c3+eps)**iexp
!        coefficients gamma en maillage irregulier
                g1p=cmuj2(m1)*cvj(m1)/(cmuj1(m)*cvj(m)+cmuj2(m)*cvj(m)+cmuj2(m1)*cvj(m1))
                g2p=1.-g1p
                g2m=cmuj1(m1)*cvj(m1)/(cmuj1(m1)*cvj(m1)+cmuj2(m1)*cvj(m1)+cmuj2(m+2*ninc)*cvj(m+2*ninc))
                g1m=1.-g2m
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
!
                ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
!
                ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
!
                ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                     +0.5*(1.-sign(1.,v4))*(g1m/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                     +0.5*(1.-sign(1.,v4))*(g2m/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
!
                ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                     +0.5*(1.-sign(1.,v5))*(g1m/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                     +0.5*(1.-sign(1.,v5))*(g2m/beta52)
                sw=ww15+ww25
                w15=ww15/sw
                w25=ww25/sw
!        calcul des flux convectifs projetes
                fc1=w11*f11+w21*f12
                fc2=w12*f21+w22*f22
                fc3=w13*f31+w23*f32
                fc4=w14*f41+w24*f42
                fc5=w15*f51+w25*f52
!        produit avec matrice P pour retour dans l'espace physique
                f1=fc1*p11+fc2*p12+fc3*p13+fc4*p14+fc5*p15
                f2=fc1*p21+fc2*p22+fc3*p23+fc4*p24+fc5*p25
                f3=fc1*p31+fc2*p32+fc3*p33+fc4*p34+fc5*p35
                f5=fc1*p51+fc2*p52+fc3*p53+fc4*p54+fc5*p55
!-----------------------------------------------------------------
!        produit de Q avec les flux Euler aux points (i-1) a (i+2)
                q1f1m=q11*v(n-ninc  ,3)+q12*fxy(m-ninc)  +q13*fyy(m-ninc) &
                     +q14*fyz(m-ninc)  +q15*fey(m-ninc)
                q1f  =q11*v(n       ,3)+q12*fxy(m)       +q13*fyy(m) &
                     +q14*fyz(m)       +q15*fey(m)
                q1f1p=q11*v(n+ninc  ,3)+q12*fxy(m+ninc)  +q13*fyy(m+ninc) &
                     +q14*fyz(m+ninc)  +q15*fey(m+ninc)
                q1f2p=q11*v(n+2*ninc,3)+q12*fxy(m+2*ninc)+q13*fyy(m+2*ninc) &
                     +q14*fyz(m+2*ninc)+q15*fey(m+2*ninc)
!
                q2f1m=q21*v(n-ninc  ,3)+q22*fxy(m-ninc)  +q23*fyy(m-ninc) &
                     +q24*fyz(m-ninc)  +q25*fey(m-ninc)
                q2f  =q21*v(n       ,3)+q22*fxy(m)       +q23*fyy(m) &
                     +q24*fyz(m)       +q25*fey(m)
                q2f1p=q21*v(n+ninc  ,3)+q22*fxy(m+ninc)  +q23*fyy(m+ninc) &
                     +q24*fyz(m+ninc)  +q25*fey(m+ninc)
                q2f2p=q21*v(n+2*ninc,3)+q22*fxy(m+2*ninc)+q23*fyy(m+2*ninc) &
                     +q24*fyz(m+2*ninc)+q25*fey(m+2*ninc)
!
                q3f1m=q31*v(n-ninc  ,3)+q32*fxy(m-ninc)  +q33*fyy(m-ninc) &
                     +q34*fyz(m-ninc)  +q35*fey(m-ninc)
                q3f  =q31*v(n       ,3)+q32*fxy(m)       +q33*fyy(m) &
                     +q34*fyz(m)       +q35*fey(m)
                q3f1p=q31*v(n+ninc  ,3)+q32*fxy(m+ninc)  +q33*fyy(m+ninc) &
                     +q34*fyz(m+ninc)  +q35*fey(m+ninc)
                q3f2p=q31*v(n+2*ninc,3)+q32*fxy(m+2*ninc)+q33*fyy(m+2*ninc) &
                     +q34*fyz(m+2*ninc)+q35*fey(m+2*ninc)
!
                q4f1m=q41*v(n-ninc  ,3)+q42*fxy(m-ninc)  +q43*fyy(m-ninc) &
                     +q44*fyz(m-ninc)  +q45*fey(m-ninc)
                q4f  =q41*v(n       ,3)+q42*fxy(m)       +q43*fyy(m) &
                     +q44*fyz(m)       +q45*fey(m)
                q4f1p=q41*v(n+ninc  ,3)+q42*fxy(m+ninc)  +q43*fyy(m+ninc) &
                     +q44*fyz(m+ninc)  +q45*fey(m+ninc)
                q4f2p=q41*v(n+2*ninc,3)+q42*fxy(m+2*ninc)+q43*fyy(m+2*ninc) &
                     +q44*fyz(m+2*ninc)+q45*fey(m+2*ninc)
!
                q5f1m=q51*v(n-ninc  ,3)+q52*fxy(m-ninc)  +q53*fyy(m-ninc) &
                     +q54*fyz(m-ninc)  +q55*fey(m-ninc)
                q5f  =q51*v(n       ,3)+q52*fxy(m)       +q53*fyy(m) &
                     +q54*fyz(m)       +q55*fey(m)
                q5f1p=q51*v(n+ninc  ,3)+q52*fxy(m+ninc)  +q53*fyy(m+ninc) &
                     +q54*fyz(m+ninc)  +q55*fey(m+ninc)
                q5f2p=q51*v(n+2*ninc,3)+q52*fxy(m+2*ninc)+q53*fyy(m+2*ninc) &
                     +q54*fyz(m+2*ninc)+q55*fey(m+2*ninc)
!        calcul des flux d'ordre 2 sur les 2 stencils
                g11=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q1f  *c00 +q1f1p*c01)
                g12=0.5*(1.+sign(1.,v1))*(q1f  *c00 +q1f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q1f1p*c11 +q1f2p*c10)
!
                g21=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q2f  *c00 +q2f1p*c01)
                g22=0.5*(1.+sign(1.,v1))*(q2f  *c00 +q2f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q2f1p*c11 +q2f2p*c10)
!
                g31=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q3f  *c00 +q3f1p*c01)
                g32=0.5*(1.+sign(1.,v1))*(q3f  *c00 +q3f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q3f1p*c11 +q3f2p*c10)
!
                g41=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11) &
                     +0.5*(1.-sign(1.,v4))*(q4f  *c00 +q4f1p*c01)
                g42=0.5*(1.+sign(1.,v4))*(q4f  *c00 +q4f1p*c01) &
                     +0.5*(1.-sign(1.,v4))*(q4f1p*c11 +q4f2p*c10)
!
                g51=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11) &
                     +0.5*(1.-sign(1.,v5))*(q5f  *c00 +q5f1p*c01)
                g52=0.5*(1.+sign(1.,v5))*(q5f  *c00 +q5f1p*c01) &
                     +0.5*(1.-sign(1.,v5))*(q5f1p*c11 +q5f2p*c10)
!        calcul des senseurs beta (au carre)
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*d1c1 &
                     +0.5*(1.-sign(1.,v1))*d1c2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*d1c2 &
                     +0.5*(1.-sign(1.,v1))*d1c3+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*d2c1 &
                     +0.5*(1.-sign(1.,v1))*d2c2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*d2c2 &
                     +0.5*(1.-sign(1.,v1))*d2c3+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*d3c1 &
                     +0.5*(1.-sign(1.,v1))*d3c2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*d3c2 &
                     +0.5*(1.-sign(1.,v1))*d3c3+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*d4c1 &
                     +0.5*(1.-sign(1.,v4))*d4c2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*d4c2 &
                     +0.5*(1.-sign(1.,v4))*d4c3+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*d5c1 &
                     +0.5*(1.-sign(1.,v5))*d5c2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*d5c2 &
                     +0.5*(1.-sign(1.,v5))*d5c3+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
!
                ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
!
                ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
!
                ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                     +0.5*(1.-sign(1.,v4))*(g1m/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                     +0.5*(1.-sign(1.,v4))*(g2m/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
!
                ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                     +0.5*(1.-sign(1.,v5))*(g1m/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                     +0.5*(1.-sign(1.,v5))*(g2m/beta52)
                sw=ww15+ww25
                w15=ww15/sw
                w25=ww25/sw
!        calcul des flux convectifs projetes
                gc1=w11*g11+w21*g12
                gc2=w12*g21+w22*g22
                gc3=w13*g31+w23*g32
                gc4=w14*g41+w24*g42
                gc5=w15*g51+w25*g52
!        produit avec matrice P pour retour dans l'espace physique
                g1=gc1*p11+gc2*p12+gc3*p13+gc4*p14+gc5*p15
                g2=gc1*p21+gc2*p22+gc3*p23+gc4*p24+gc5*p25
                g3=gc1*p31+gc2*p32+gc3*p33+gc4*p34+gc5*p35
                g5=gc1*p51+gc2*p52+gc3*p53+gc4*p54+gc5*p55
!        calcul du flux numerique et bilan de flux
                dg1=f1*sn(m1,kdir,1)+g1*sn(m1,kdir,2)
                dg2=f2*sn(m1,kdir,1)+g2*sn(m1,kdir,2)
                dg3=f3*sn(m1,kdir,1)+g3*sn(m1,kdir,2)
                dg5=f5*sn(m1,kdir,1)+g5*sn(m1,kdir,2)
!        calcul des flux visqueux (multiplies par -2)
                gv2=(cmuj2(m1)*toxx(n)+cmuj1(m1)*toxx(n1))*sn(m1,kdir,1) &
                     +(cmuj2(m1)*toxy(n)+cmuj1(m1)*toxy(n1))*sn(m1,kdir,2)
                gv3=(cmuj2(m1)*toxy(n)+cmuj1(m1)*toxy(n1))*sn(m1,kdir,1) &
                     +(cmuj2(m1)*toyy(n)+cmuj1(m1)*toyy(n1))*sn(m1,kdir,2)
                gv5=(cmuj2(m1)*(toxx(n )*ul+toxy(n )*vl+qcx(n )) &
                     +cmuj1(m1)*(toxx(n1)*ur+toxy(n1)*vr+qcx(n1)))*sn(m1,kdir,1) &
                     +(cmuj2(m1)*(toxy(n )*ul+toyy(n )*vl+qcy(n )) &
                     +cmuj1(m1)*(toxy(n1)*ur+toyy(n1)*vr+qcy(n1)))*sn(m1,kdir,2)
                u(n1,1)=u(n1,1)-dg1
                u(n1,2)=u(n1,2)-dg2+0.5*gv2
                u(n1,3)=u(n1,3)-dg3+0.5*gv3
                u(n1,5)=u(n1,5)-dg5+0.5*gv5
                u(n,1)=u(n,1)+dg1
                u(n,2)=u(n,2)+dg2-0.5*gv2
                u(n,3)=u(n,3)+dg3-0.5*gv3
                u(n,5)=u(n,5)+dg5-0.5*gv5
             enddo
          enddo
       enddo

    elseif(imap.eq.1) then
!
!-----direction i-----------------------------------------
!
       kdir=1
       ninc=nci
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m2,j,k)
!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                m1=m+ninc
                n1=n+ninc
!        vecteur normal unitaire a la face consideree (face i+1/2)
                cnds=sqrt(sn(m1,kdir,1)*sn(m1,kdir,1)+ &
                     sn(m1,kdir,2)*sn(m1,kdir,2)+ &
                     sn(m1,kdir,3)*sn(m1,kdir,3))
                nx=sn(m1,kdir,1)/cnds
                ny=sn(m1,kdir,2)/cnds
                nz=sn(m1,kdir,3)/cnds
!        calcul des etats gauche et droit
                ul=v(n,2)/v(n,1)
                vl=v(n,3)/v(n,1)
                wl=v(n,4)/v(n,1)
                ur=v(n1,2)/v(n1,1)
                vr=v(n1,3)/v(n1,1)
                wr=v(n1,4)/v(n1,1)
                al=sqrt(gam*ps(n )/v(n,1))
                ar=sqrt(gam*ps(n1)/v(n1,1))
                hl=al*al/gam1+0.5*(ul**2+vl**2+wl**2)
                hr=ar*ar/gam1+0.5*(ur**2+vr**2+wr**2)
!        calcul de etat moyen de Roe
                gd=sqrt(v(n1,1)/v(n,1))
                gd1=1./(1.+gd)
                gd2=gd*gd1
                rhom=sqrt(v(n,1)*v(n1,1))
                rhomi=1./rhom
                um=gd1*ul+gd2*ur
                vm=gd1*vl+gd2*vr
                wm=gd1*wl+gd2*wr
                hm=gd1*hl+gd2*hr
                vitm2=0.5*(um**2+vm**2+wm**2)
                am=sqrt(abs(gam1*(hm-vitm2)))
                am2i=1./(am*am)
                vn=um*nx+vm*ny+wm*nz
                rhoiam=rhom/am
                rhoami=am2i/rhoiam
!        valeurs propres de la matrice jacobienne des flux
                v1=vn
                v4=vn+am
                v5=vn-am
!        calcul des matrices de passage a gauche Q et a droite P
                q11=(1.-gam1*vitm2*am2i)*nx-(vm*nz-wm*ny)*rhomi
                q12=gam1*um*nx*am2i
                q13=gam1*vm*nx*am2i+nz*rhomi
                q14=gam1*wm*nx*am2i-ny*rhomi
                q15=-gam1*nx*am2i
                q21=(1.-gam1*vitm2*am2i)*ny-(wm*nx-um*nz)*rhomi
                q22=gam1*um*ny*am2i-nz*rhomi
                q23=gam1*vm*ny*am2i
                q24=gam1*wm*ny*am2i+nx*rhomi
                q25=-gam1*ny*am2i
                q31=(1.-gam1*vitm2*am2i)*nz-(um*ny-vm*nx)*rhomi
                q32=gam1*um*nz*am2i+ny*rhomi
                q33=gam1*vm*nz*am2i-nx*rhomi
                q34=gam1*wm*nz*am2i
                q35=-gam1*nz*am2i
                q41=gam1*vitm2*rhoami-vn*rhomi
                q42=nx*rhomi-gam1*um*rhoami
                q43=ny*rhomi-gam1*vm*rhoami
                q44=nz*rhomi-gam1*wm*rhoami
                q45=gam1*rhoami
                q51=gam1*vitm2*rhoami+vn*rhomi
                q52=-nx*rhomi-gam1*um*rhoami
                q53=-ny*rhomi-gam1*vm*rhoami
                q54=-nz*rhomi-gam1*wm*rhoami
                q55=gam1*rhoami
!
                p11=nx
                p12=ny
                p13=nz
                p14=0.5*rhoiam
                p15=0.5*rhoiam
                p21=um*nx
                p22=um*ny-rhom*nz
                p23=um*nz+rhom*ny
                p24=0.5*rhoiam*(um+nx*am)
                p25=0.5*rhoiam*(um-nx*am)
                p31=vm*nx+rhom*nz
                p32=vm*ny
                p33=vm*nz-rhom*nx
                p34=0.5*rhoiam*(vm+ny*am)
                p35=0.5*rhoiam*(vm-ny*am)
                p41=wm*nx-rhom*ny
                p42=wm*ny+rhom*nx
                p43=wm*nz
                p44=0.5*rhoiam*(wm+nz*am)
                p45=0.5*rhoiam*(wm-nz*am)
                p51=vitm2*nx+rhom*(vm*nz-wm*ny)
                p52=vitm2*ny+rhom*(wm*nx-um*nz)
                p53=vitm2*nz+rhom*(um*ny-vm*nx)
                p54=0.5*rhoiam*(hm+am*vn)
                p55=0.5*rhoiam*(hm-am*vn)
!        produit de Q avec les flux Euler aux points (i-1) a (i+2)
                q1f1m=q11*v(n-ninc  ,2)+q12*fxx(m-ninc)  +q13*fxy(m-ninc) &
                     +q14*fxz(m-ninc)  +q15*fex(m-ninc)
                q1f  =q11*v(n       ,2)+q12*fxx(m)       +q13*fxy(m) &
                     +q14*fxz(m)       +q15*fex(m)
                q1f1p=q11*v(n+ninc  ,2)+q12*fxx(m+ninc)  +q13*fxy(m+ninc) &
                     +q14*fxz(m+ninc)  +q15*fex(m+ninc)
                q1f2p=q11*v(n+2*ninc,2)+q12*fxx(m+2*ninc)+q13*fxy(m+2*ninc) &
                     +q14*fxz(m+2*ninc)+q15*fex(m+2*ninc)
!
                q2f1m=q21*v(n-ninc  ,2)+q22*fxx(m-ninc)  +q23*fxy(m-ninc) &
                     +q24*fxz(m-ninc)  +q25*fex(m-ninc)
                q2f  =q21*v(n       ,2)+q22*fxx(m)       +q23*fxy(m) &
                     +q24*fxz(m)       +q25*fex(m)
                q2f1p=q21*v(n+ninc  ,2)+q22*fxx(m+ninc)  +q23*fxy(m+ninc) &
                     +q24*fxz(m+ninc)  +q25*fex(m+ninc)
                q2f2p=q21*v(n+2*ninc,2)+q22*fxx(m+2*ninc)+q23*fxy(m+2*ninc) &
                     +q24*fxz(m+2*ninc)+q25*fex(m+2*ninc)
!
                q3f1m=q31*v(n-ninc  ,2)+q32*fxx(m-ninc)  +q33*fxy(m-ninc) &
                     +q34*fxz(m-ninc)  +q35*fex(m-ninc)
                q3f  =q31*v(n       ,2)+q32*fxx(m)       +q33*fxy(m) &
                     +q34*fxz(m)       +q35*fex(m)
                q3f1p=q31*v(n+ninc  ,2)+q32*fxx(m+ninc)  +q33*fxy(m+ninc) &
                     +q34*fxz(m+ninc)  +q35*fex(m+ninc)
                q3f2p=q31*v(n+2*ninc,2)+q32*fxx(m+2*ninc)+q33*fxy(m+2*ninc) &
                     +q34*fxz(m+2*ninc)+q35*fex(m+2*ninc)
!
                q4f1m=q41*v(n-ninc  ,2)+q42*fxx(m-ninc)  +q43*fxy(m-ninc) &
                     +q44*fxz(m-ninc)  +q45*fex(m-ninc)
                q4f  =q41*v(n       ,2)+q42*fxx(m)       +q43*fxy(m) &
                     +q44*fxz(m)       +q45*fex(m)
                q4f1p=q41*v(n+ninc  ,2)+q42*fxx(m+ninc)  +q43*fxy(m+ninc) &
                     +q44*fxz(m+ninc)  +q45*fex(m+ninc)
                q4f2p=q41*v(n+2*ninc,2)+q42*fxx(m+2*ninc)+q43*fxy(m+2*ninc) &
                     +q44*fxz(m+2*ninc)+q45*fex(m+2*ninc)
!
                q5f1m=q51*v(n-ninc  ,2)+q52*fxx(m-ninc)  +q53*fxy(m-ninc) &
                     +q54*fxz(m-ninc)  +q55*fex(m-ninc)
                q5f  =q51*v(n       ,2)+q52*fxx(m)       +q53*fxy(m) &
                     +q54*fxz(m)       +q55*fex(m)
                q5f1p=q51*v(n+ninc  ,2)+q52*fxx(m+ninc)  +q53*fxy(m+ninc) &
                     +q54*fxz(m+ninc)  +q55*fex(m+ninc)
                q5f2p=q51*v(n+2*ninc,2)+q52*fxx(m+2*ninc)+q53*fxy(m+2*ninc) &
                     +q54*fxz(m+2*ninc)+q55*fex(m+2*ninc)
!        coefficients Crj en maillage irregulier
                c00=cmui2(m1)*0.5
                c01=cmui1(m1)*0.5
                c10=-cmui2(m)*0.5
                c11=1.-c10
                c21=-cmui1(m+2*ninc)*0.5
                c20=1.-c21
!        calcul des flux d'ordre 2 sur les 2 stencils
                f11=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q1f  *c00 +q1f1p*c01)
                f12=0.5*(1.+sign(1.,v1))*(q1f  *c00 +q1f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q1f1p*c11 +q1f2p*c10)
!
                f21=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q2f  *c00 +q2f1p*c01)
                f22=0.5*(1.+sign(1.,v1))*(q2f  *c00 +q2f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q2f1p*c11 +q2f2p*c10)
!
                f31=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q3f  *c00 +q3f1p*c01)
                f32=0.5*(1.+sign(1.,v1))*(q3f  *c00 +q3f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q3f1p*c11 +q3f2p*c10)
!
                f41=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11) &
                     +0.5*(1.-sign(1.,v4))*(q4f  *c00 +q4f1p*c01)
                f42=0.5*(1.+sign(1.,v4))*(q4f  *c00 +q4f1p*c01) &
                     +0.5*(1.-sign(1.,v4))*(q4f1p*c11 +q4f2p*c10)
!
                f51=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11) &
                     +0.5*(1.-sign(1.,v5))*(q5f  *c00 +q5f1p*c01)
                f52=0.5*(1.+sign(1.,v5))*(q5f  *c00 +q5f1p*c01) &
                     +0.5*(1.-sign(1.,v5))*(q5f1p*c11 +q5f2p*c10)
!        calcul des senseurs beta (au carre)
!         d1c1=(q1f-q1f1m)**2    !maillage regulier
!         d1c2=(q1f1p-q1f)**2
!         d1c3=(q1f2p-q1f1p)**2
                d1c1=4.*((q1f-q1f1m)*c10)**2
                d1c2=4.*((q1f1p-q1f)*c01)**2
                d1c3=4.*((q1f2p-q1f1p)*c01*c21/c00)**2
                d2c1=4.*((q2f-q2f1m)*c10)**2
                d2c2=4.*((q2f1p-q2f)*c01)**2
                d2c3=4.*((q2f2p-q2f1p)*c01*c21/c00)**2
                d3c1=4.*((q3f-q3f1m)*c10)**2
                d3c2=4.*((q3f1p-q3f)*c01)**2
                d3c3=4.*((q3f2p-q3f1p)*c01*c21/c00)**2
                d4c1=4.*((q4f-q4f1m)*c10)**2
                d4c2=4.*((q4f1p-q4f)*c01)**2
                d4c3=4.*((q4f2p-q4f1p)*c01*c21/c00)**2
                d5c1=4.*((q5f-q5f1m)*c10)**2
                d5c2=4.*((q5f1p-q5f)*c01)**2
                d5c3=4.*((q5f2p-q5f1p)*c01*c21/c00)**2
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*d1c1 &
                     +0.5*(1.-sign(1.,v1))*d1c2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*d1c2 &
                     +0.5*(1.-sign(1.,v1))*d1c3+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*d2c1 &
                     +0.5*(1.-sign(1.,v1))*d2c2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*d2c2 &
                     +0.5*(1.-sign(1.,v1))*d2c3+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*d3c1 &
                     +0.5*(1.-sign(1.,v1))*d3c2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*d3c2 &
                     +0.5*(1.-sign(1.,v1))*d3c3+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*d4c1 &
                     +0.5*(1.-sign(1.,v4))*d4c2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*d4c2 &
                     +0.5*(1.-sign(1.,v4))*d4c3+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*d5c1 &
                     +0.5*(1.-sign(1.,v5))*d5c2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*d5c2 &
                     +0.5*(1.-sign(1.,v5))*d5c3+eps)**iexp
!        coefficients gamma en maillage irregulier
                g1p=cmui2(m1)*cvi(m1)/(cmui1(m)*cvi(m)+cmui2(m)*cvi(m)+cmui2(m1)*cvi(m1))
                g2p=1.-g1p
                g2m=cmui1(m1)*cvi(m1)/(cmui1(m1)*cvi(m1)+cmui2(m1)*cvi(m1)+cmui2(m+2*ninc)*cvi(m+2*ninc))
                g1m=1.-g2m
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
                ww11m=w11*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w11+w11**2)/(g1p**2+w11*(1.-2.*g1p)) &
                     +w11*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w11+w11**2)/(g1m**2+w11*(1.-2.*g1m))
                ww21m=w21*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w21+w21**2)/(g2p**2+w21*(1.-2.*g2p)) &
                     +w21*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w21+w21**2)/(g2m**2+w21*(1.-2.*g2m))
                swm=ww11m+ww21m
                w11=ww11m/swm
                w21=ww21m/swm
!
                ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
                ww12m=w12*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w12+w12**2)/(g1p**2+w12*(1.-2.*g1p)) &
                     +w12*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w12+w12**2)/(g1m**2+w12*(1.-2.*g1m))
                ww22m=w22*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w22+w22**2)/(g2p**2+w22*(1.-2.*g2p)) &
                     +w22*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w22+w22**2)/(g2m**2+w22*(1.-2.*g2m))
                swm=ww12m+ww22m
                w12=ww12m/swm
                w22=ww22m/swm
!
                ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
                ww13m=w13*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w13+w13**2)/(g1p**2+w13*(1.-2.*g1p)) &
                     +w13*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w13+w13**2)/(g1m**2+w13*(1.-2.*g1m))
                ww23m=w23*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w23+w23**2)/(g2p**2+w23*(1.-2.*g2p)) &
                     +w23*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w23+w23**2)/(g2m**2+w23*(1.-2.*g2m))
                swm=ww13m+ww23m
                w13=ww13m/swm
                w23=ww23m/swm
!
                ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                     +0.5*(1.-sign(1.,v4))*(g1m/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                     +0.5*(1.-sign(1.,v4))*(g2m/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
                ww14m=w14*0.5*(1.+sign(1.,v4))*(g1p+g1p**2-3.*g1p*w14+w14**2)/(g1p**2+w14*(1.-2.*g1p)) &
                     +w14*0.5*(1.-sign(1.,v4))*(g1m+g1m**2-3.*g1m*w14+w14**2)/(g1m**2+w14*(1.-2.*g1m))
                ww24m=w24*0.5*(1.+sign(1.,v4))*(g2p+g2p**2-3.*g2p*w24+w24**2)/(g2p**2+w24*(1.-2.*g2p)) &
                     +w24*0.5*(1.-sign(1.,v4))*(g2m+g2m**2-3.*g2m*w24+w24**2)/(g2m**2+w24*(1.-2.*g2m))
                swm=ww14m+ww24m
                w14=ww14m/swm
                w24=ww24m/swm
!
                ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                     +0.5*(1.-sign(1.,v5))*(g1m/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                     +0.5*(1.-sign(1.,v5))*(g2m/beta52)
                sw=ww15+ww25
                w15=ww15/sw
                w25=ww25/sw
                ww15m=w15*0.5*(1.+sign(1.,v5))*(g1p+g1p**2-3.*g1p*w15+w15**2)/(g1p**2+w15*(1.-2.*g1p)) &
                     +w15*0.5*(1.-sign(1.,v5))*(g1m+g1m**2-3.*g1m*w15+w15**2)/(g1m**2+w15*(1.-2.*g1m))
                ww25m=w25*0.5*(1.+sign(1.,v5))*(g2p+g2p**2-3.*g2p*w25+w25**2)/(g2p**2+w25*(1.-2.*g2p)) &
                     +w25*0.5*(1.-sign(1.,v5))*(g2m+g2m**2-3.*g2m*w25+w25**2)/(g2m**2+w25*(1.-2.*g2m))
                swm=ww15m+ww25m
                w15=ww15m/swm
                w25=ww25m/swm
!        calcul des flux convectifs projetes
                fc1=w11*f11+w21*f12
                fc2=w12*f21+w22*f22
                fc3=w13*f31+w23*f32
                fc4=w14*f41+w24*f42
                fc5=w15*f51+w25*f52
!        produit avec matrice P pour retour dans l'espace physique
                f1=fc1*p11+fc2*p12+fc3*p13+fc4*p14+fc5*p15
                f2=fc1*p21+fc2*p22+fc3*p23+fc4*p24+fc5*p25
                f3=fc1*p31+fc2*p32+fc3*p33+fc4*p34+fc5*p35
                f5=fc1*p51+fc2*p52+fc3*p53+fc4*p54+fc5*p55
!-------------------------------------------------------------------
!        produit de Q avec les flux Euler aux points (i-1) a (i+2)
                q1f1m=q11*v(n-ninc  ,3)+q12*fxy(m-ninc)  +q13*fyy(m-ninc) &
                     +q14*fyz(m-ninc)  +q15*fey(m-ninc)
                q1f  =q11*v(n       ,3)+q12*fxy(m)       +q13*fyy(m) &
                     +q14*fyz(m)       +q15*fey(m)
                q1f1p=q11*v(n+ninc  ,3)+q12*fxy(m+ninc)  +q13*fyy(m+ninc) &
                     +q14*fyz(m+ninc)  +q15*fey(m+ninc)
                q1f2p=q11*v(n+2*ninc,3)+q12*fxy(m+2*ninc)+q13*fyy(m+2*ninc) &
                     +q14*fyz(m+2*ninc)+q15*fey(m+2*ninc)
!
                q2f1m=q21*v(n-ninc  ,3)+q22*fxy(m-ninc)  +q23*fyy(m-ninc) &
                     +q24*fyz(m-ninc)  +q25*fey(m-ninc)
                q2f  =q21*v(n       ,3)+q22*fxy(m)       +q23*fyy(m) &
                     +q24*fyz(m)       +q25*fey(m)
                q2f1p=q21*v(n+ninc  ,3)+q22*fxy(m+ninc)  +q23*fyy(m+ninc) &
                     +q24*fyz(m+ninc)  +q25*fey(m+ninc)
                q2f2p=q21*v(n+2*ninc,3)+q22*fxy(m+2*ninc)+q23*fyy(m+2*ninc) &
                     +q24*fyz(m+2*ninc)+q25*fey(m+2*ninc)
!
                q3f1m=q31*v(n-ninc  ,3)+q32*fxy(m-ninc)  +q33*fyy(m-ninc) &
                     +q34*fyz(m-ninc)  +q35*fey(m-ninc)
                q3f  =q31*v(n       ,3)+q32*fxy(m)       +q33*fyy(m) &
                     +q34*fyz(m)       +q35*fey(m)
                q3f1p=q31*v(n+ninc  ,3)+q32*fxy(m+ninc)  +q33*fyy(m+ninc) &
                     +q34*fyz(m+ninc)  +q35*fey(m+ninc)
                q3f2p=q31*v(n+2*ninc,3)+q32*fxy(m+2*ninc)+q33*fyy(m+2*ninc) &
                     +q34*fyz(m+2*ninc)+q35*fey(m+2*ninc)
!
                q4f1m=q41*v(n-ninc  ,3)+q42*fxy(m-ninc)  +q43*fyy(m-ninc) &
                     +q44*fyz(m-ninc)  +q45*fey(m-ninc)
                q4f  =q41*v(n       ,3)+q42*fxy(m)       +q43*fyy(m) &
                     +q44*fyz(m)       +q45*fey(m)
                q4f1p=q41*v(n+ninc  ,3)+q42*fxy(m+ninc)  +q43*fyy(m+ninc) &
                     +q44*fyz(m+ninc)  +q45*fey(m+ninc)
                q4f2p=q41*v(n+2*ninc,3)+q42*fxy(m+2*ninc)+q43*fyy(m+2*ninc) &
                     +q44*fyz(m+2*ninc)+q45*fey(m+2*ninc)
!
                q5f1m=q51*v(n-ninc  ,3)+q52*fxy(m-ninc)  +q53*fyy(m-ninc) &
                     +q54*fyz(m-ninc)  +q55*fey(m-ninc)
                q5f  =q51*v(n       ,3)+q52*fxy(m)       +q53*fyy(m) &
                     +q54*fyz(m)       +q55*fey(m)
                q5f1p=q51*v(n+ninc  ,3)+q52*fxy(m+ninc)  +q53*fyy(m+ninc) &
                     +q54*fyz(m+ninc)  +q55*fey(m+ninc)
                q5f2p=q51*v(n+2*ninc,3)+q52*fxy(m+2*ninc)+q53*fyy(m+2*ninc) &
                     +q54*fyz(m+2*ninc)+q55*fey(m+2*ninc)
!        calcul des flux d'ordre 2 sur les 2 stencils
                g11=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q1f  *c00 +q1f1p*c01)
                g12=0.5*(1.+sign(1.,v1))*(q1f  *c00 +q1f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q1f1p*c11 +q1f2p*c10)
!
                g21=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q2f  *c00 +q2f1p*c01)
                g22=0.5*(1.+sign(1.,v1))*(q2f  *c00 +q2f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q2f1p*c11 +q2f2p*c10)
!
                g31=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q3f  *c00 +q3f1p*c01)
                g32=0.5*(1.+sign(1.,v1))*(q3f  *c00 +q3f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q3f1p*c11 +q3f2p*c10)
!
                g41=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11) &
                     +0.5*(1.-sign(1.,v4))*(q4f  *c00 +q4f1p*c01)
                g42=0.5*(1.+sign(1.,v4))*(q4f  *c00 +q4f1p*c01) &
                     +0.5*(1.-sign(1.,v4))*(q4f1p*c11 +q4f2p*c10)
!
                g51=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11) &
                     +0.5*(1.-sign(1.,v5))*(q5f  *c00 +q5f1p*c01)
                g52=0.5*(1.+sign(1.,v5))*(q5f  *c00 +q5f1p*c01) &
                     +0.5*(1.-sign(1.,v5))*(q5f1p*c11 +q5f2p*c10)
!        calcul des senseurs beta (au carre)
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*d1c1 &
                     +0.5*(1.-sign(1.,v1))*d1c2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*d1c2 &
                     +0.5*(1.-sign(1.,v1))*d1c3+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*d2c1 &
                     +0.5*(1.-sign(1.,v1))*d2c2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*d2c2 &
                     +0.5*(1.-sign(1.,v1))*d2c3+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*d3c1 &
                     +0.5*(1.-sign(1.,v1))*d3c2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*d3c2 &
                     +0.5*(1.-sign(1.,v1))*d3c3+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*d4c1 &
                     +0.5*(1.-sign(1.,v4))*d4c2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*d4c2 &
                     +0.5*(1.-sign(1.,v4))*d4c3+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*d5c1 &
                     +0.5*(1.-sign(1.,v5))*d5c2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*d5c2 &
                     +0.5*(1.-sign(1.,v5))*d5c3+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
                ww11m=w11*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w11+w11**2)/(g1p**2+w11*(1.-2.*g1p)) &
                     +w11*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w11+w11**2)/(g1m**2+w11*(1.-2.*g1m))
                ww21m=w21*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w21+w21**2)/(g2p**2+w21*(1.-2.*g2p)) &
                     +w21*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w21+w21**2)/(g2m**2+w21*(1.-2.*g2m))
                swm=ww11m+ww21m
                w11=ww11m/swm
                w21=ww21m/swm
!
                ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
                ww12m=w12*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w12+w12**2)/(g1p**2+w12*(1.-2.*g1p)) &
                     +w12*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w12+w12**2)/(g1m**2+w12*(1.-2.*g1m))
                ww22m=w22*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w22+w22**2)/(g2p**2+w22*(1.-2.*g2p)) &
                     +w22*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w22+w22**2)/(g2m**2+w22*(1.-2.*g2m))
                swm=ww12m+ww22m
                w12=ww12m/swm
                w22=ww22m/swm
!
                ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
                ww13m=w13*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w13+w13**2)/(g1p**2+w13*(1.-2.*g1p)) &
                     +w13*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w13+w13**2)/(g1m**2+w13*(1.-2.*g1m))
                ww23m=w23*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w23+w23**2)/(g2p**2+w23*(1.-2.*g2p)) &
                     +w23*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w23+w23**2)/(g2m**2+w23*(1.-2.*g2m))
                swm=ww13m+ww23m
                w13=ww13m/swm
                w23=ww23m/swm
!
                ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                     +0.5*(1.-sign(1.,v4))*(g1m/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                     +0.5*(1.-sign(1.,v4))*(g2m/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
                ww14m=w14*0.5*(1.+sign(1.,v4))*(g1p+g1p**2-3.*g1p*w14+w14**2)/(g1p**2+w14*(1.-2.*g1p)) &
                     +w14*0.5*(1.-sign(1.,v4))*(g1m+g1m**2-3.*g1m*w14+w14**2)/(g1m**2+w14*(1.-2.*g1m))
                ww24m=w24*0.5*(1.+sign(1.,v4))*(g2p+g2p**2-3.*g2p*w24+w24**2)/(g2p**2+w24*(1.-2.*g2p)) &
                     +w24*0.5*(1.-sign(1.,v4))*(g2m+g2m**2-3.*g2m*w24+w24**2)/(g2m**2+w24*(1.-2.*g2m))
                swm=ww14m+ww24m
                w14=ww14m/swm
                w24=ww24m/swm
!
                ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                     +0.5*(1.-sign(1.,v5))*(g1m/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                     +0.5*(1.-sign(1.,v5))*(g2m/beta52)
                sw=ww15+ww25
                w15=ww15/sw
                w25=ww25/sw
                ww15m=w15*0.5*(1.+sign(1.,v5))*(g1p+g1p**2-3.*g1p*w15+w15**2)/(g1p**2+w15*(1.-2.*g1p)) &
                     +w15*0.5*(1.-sign(1.,v5))*(g1m+g1m**2-3.*g1m*w15+w15**2)/(g1m**2+w15*(1.-2.*g1m))
                ww25m=w25*0.5*(1.+sign(1.,v5))*(g2p+g2p**2-3.*g2p*w25+w25**2)/(g2p**2+w25*(1.-2.*g2p)) &
                     +w25*0.5*(1.-sign(1.,v5))*(g2m+g2m**2-3.*g2m*w25+w25**2)/(g2m**2+w25*(1.-2.*g2m))
                swm=ww15m+ww25m
                w15=ww15m/swm
                w25=ww25m/swm
!        calcul des flux convectifs projetes
                gc1=w11*g11+w21*g12
                gc2=w12*g21+w22*g22
                gc3=w13*g31+w23*g32
                gc4=w14*g41+w24*g42
                gc5=w15*g51+w25*g52
!        produit avec matrice P pour retour dans l'espace physique
                g1=gc1*p11+gc2*p12+gc3*p13+gc4*p14+gc5*p15
                g2=gc1*p21+gc2*p22+gc3*p23+gc4*p24+gc5*p25
                g3=gc1*p31+gc2*p32+gc3*p33+gc4*p34+gc5*p35
                g5=gc1*p51+gc2*p52+gc3*p53+gc4*p54+gc5*p55
!        calcul du flux numerique et bilan de flux
                df1=f1*sn(m1,kdir,1)+g1*sn(m1,kdir,2)
                df2=f2*sn(m1,kdir,1)+g2*sn(m1,kdir,2)
                df3=f3*sn(m1,kdir,1)+g3*sn(m1,kdir,2)
                df5=f5*sn(m1,kdir,1)+g5*sn(m1,kdir,2)
!        correction point sonique
!         vnl=ul*nx+vl*ny+wl*nz
!         vnr=ur*nx+vr*ny+wr*nz
!         if((vnl-al)*(vnl-ar).lt.0.) then
!          rhoel=v(n,1 )*hl-(ps(n) -pinfl)
!          rhoer=v(n1,1)*hr-(ps(n1)-pinfl)
!          rm=cnds*max(abs(vnl+al),abs(vnr+ar))
!          di5=rm*(rhoer-rhoel)
!!          dw1=v(n1,1)   -v(n,1)
!!          dw2=v(n1,1)*ur-v(n,1)*ul
!!          dw3=v(n1,1)*vr-v(n,1)*vl
!!          dw4=v(n1,1)*wr-v(n,1)*wl
!!          dw5=rhoer     -rhoel
!!          di5=(p51*v1*q11+p52*v1*q21+p53*v1*q31+p54*v4*q41+p55*v5*q51)*dw1 &
!!             +(p51*v1*q12+p52*v1*q22+p53*v1*q32+p54*v4*q42+p55*v5*q52)*dw2 &
!!             +(p51*v1*q13+p52*v1*q23+p53*v1*q33+p54*v4*q43+p55*v5*q53)*dw3 &
!!             +(p51*v1*q14+p52*v1*q24+p53*v1*q34+p54*v4*q44+p55*v5*q54)*dw4 &
!!             +(p51*v1*q15+p52*v1*q25+p53*v1*q35+p54*v4*q45+p55*v5*q55)*dw5
!          dfex=cmui1(m)*(rhoer+ps(n1)-pinfl)*ur  &
!              +cmui2(m)*(rhoel+ps(n )-pinfl)*ul
!          dfey=cmui1(m)*(rhoer+ps(n1)-pinfl)*vr  &
!              +cmui2(m)*(rhoel+ps(n )-pinfl)*vl
!          df5=0.5*(dfex*sn(m,kdir,1)+dfey*sn(m,kdir,2)-di5)
!         endif
!        calcul des flux visqueux (multiplies par -2)
                fv2=(cmui2(m1)*toxx(n)+cmui1(m1)*toxx(n1))*sn(m1,kdir,1) &
                     +(cmui2(m1)*toxy(n)+cmui1(m1)*toxy(n1))*sn(m1,kdir,2)
                fv3=(cmui2(m1)*toxy(n)+cmui1(m1)*toxy(n1))*sn(m1,kdir,1) &
                     +(cmui2(m1)*toyy(n)+cmui1(m1)*toyy(n1))*sn(m1,kdir,2)
                fv5=(cmui2(m1)*(toxx(n )*ul+toxy(n )*vl+qcx(n )) &
                     +cmui1(m1)*(toxx(n1)*ur+toxy(n1)*vr+qcx(n1)))*sn(m1,kdir,1) &
                     +(cmui2(m1)*(toxy(n )*ul+toyy(n )*vl+qcy(n )) &
                     +cmui1(m1)*(toxy(n1)*ur+toyy(n1)*vr+qcy(n1)))*sn(m1,kdir,2)
                u(n1,1)=u(n1,1)-df1
                u(n1,2)=u(n1,2)-df2+0.5*fv2
                u(n1,3)=u(n1,3)-df3+0.5*fv3
                u(n1,5)=u(n1,5)-df5+0.5*fv5
                u(n,1)=u(n,1)+df1
                u(n,2)=u(n,2)+df2-0.5*fv2
                u(n,3)=u(n,3)+df3-0.5*fv3
                u(n,5)=u(n,5)+df5-0.5*fv5
             enddo
          enddo
       enddo
!
!------direction j----------------------------------------------
!
       kdir=2
       ninc=ncj
!
       do k=k1,k2m1
          do j=j1,j2m2
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                m1=m+ninc
                n1=n+ninc
!        vecteur normal unitaire a la face consideree (face j+1/2)
                cnds=sqrt(sn(m1,kdir,1)*sn(m1,kdir,1)+ &
                     sn(m1,kdir,2)*sn(m1,kdir,2)+ &
                     sn(m1,kdir,3)*sn(m1,kdir,3))
                nx=sn(m1,kdir,1)/cnds
                ny=sn(m1,kdir,2)/cnds
                nz=sn(m1,kdir,3)/cnds
!        calcul des etats gauche et droit
                ul=v(n,2)/v(n,1)
                vl=v(n,3)/v(n,1)
                wl=v(n,4)/v(n,1)
                ur=v(n1,2)/v(n1,1)
                vr=v(n1,3)/v(n1,1)
                wr=v(n1,4)/v(n1,1)
                al=sqrt(gam*ps(n )/v(n,1))
                ar=sqrt(gam*ps(n1)/v(n1,1))
                hl=al*al/gam1+0.5*(ul**2+vl**2+wl**2)
                hr=ar*ar/gam1+0.5*(ur**2+vr**2+wr**2)
!        calcul des etats moyens de Roe
                gd=sqrt(v(n1,1)/v(n,1))
                gd1=1./(1.+gd)
                gd2=gd*gd1
                rhom=sqrt(v(n,1)*v(n1,1))
                rhomi=1./rhom
                um=gd1*ul+gd2*ur
                vm=gd1*vl+gd2*vr
                wm=gd1*wl+gd2*wr
                hm=gd1*hl+gd2*hr
                vitm2=0.5*(um**2+vm**2+wm**2)
                am=sqrt(abs(gam1*(hm-vitm2)))
                am2i=1./(am*am)
                vn=um*nx+vm*ny+wm*nz
                rhoiam=rhom/am
                rhoami=am2i/rhoiam
!        valeurs propres
                v1=vn
                v4=vn+am
                v5=vn-am
!        calcul des matrices de passage a gauche Q et a droite P
                q11=(1.-gam1*vitm2*am2i)*nx-(vm*nz-wm*ny)*rhomi
                q12=gam1*um*nx*am2i
                q13=gam1*vm*nx*am2i+nz*rhomi
                q14=gam1*wm*nx*am2i-ny*rhomi
                q15=-gam1*nx*am2i
                q21=(1.-gam1*vitm2*am2i)*ny-(wm*nx-um*nz)*rhomi
                q22=gam1*um*ny*am2i-nz*rhomi
                q23=gam1*vm*ny*am2i
                q24=gam1*wm*ny*am2i+nx*rhomi
                q25=-gam1*ny*am2i
                q31=(1.-gam1*vitm2*am2i)*nz-(um*ny-vm*nx)*rhomi
                q32=gam1*um*nz*am2i+ny*rhomi
                q33=gam1*vm*nz*am2i-nx*rhomi
                q34=gam1*wm*nz*am2i
                q35=-gam1*nz*am2i
                q41=gam1*vitm2*rhoami-vn*rhomi
                q42=nx*rhomi-gam1*um*rhoami
                q43=ny*rhomi-gam1*vm*rhoami
                q44=nz*rhomi-gam1*wm*rhoami
                q45=gam1*rhoami
                q51=gam1*vitm2*rhoami+vn*rhomi
                q52=-nx*rhomi-gam1*um*rhoami
                q53=-ny*rhomi-gam1*vm*rhoami
                q54=-nz*rhomi-gam1*wm*rhoami
                q55=gam1*rhoami
!
                p11=nx
                p12=ny
                p13=nz
                p14=0.5*rhoiam
                p15=0.5*rhoiam
                p21=um*nx
                p22=um*ny-rhom*nz
                p23=um*nz+rhom*ny
                p24=0.5*rhoiam*(um+nx*am)
                p25=0.5*rhoiam*(um-nx*am)
                p31=vm*nx+rhom*nz
                p32=vm*ny
                p33=vm*nz-rhom*nx
                p34=0.5*rhoiam*(vm+ny*am)
                p35=0.5*rhoiam*(vm-ny*am)
                p41=wm*nx-rhom*ny
                p42=wm*ny+rhom*nx
                p43=wm*nz
                p44=0.5*rhoiam*(wm+nz*am)
                p45=0.5*rhoiam*(wm-nz*am)
                p51=vitm2*nx+rhom*(vm*nz-wm*ny)
                p52=vitm2*ny+rhom*(wm*nx-um*nz)
                p53=vitm2*nz+rhom*(um*ny-vm*nx)
                p54=0.5*rhoiam*(hm+am*vn)
                p55=0.5*rhoiam*(hm-am*vn)
!        produit de Q avec les flux Euler aux points (i-1) a (i+2)
                q1f1m=q11*v(n-ninc  ,2)+q12*fxx(m-ninc)  +q13*fxy(m-ninc) &
                     +q14*fxz(m-ninc)  +q15*fex(m-ninc)
                q1f  =q11*v(n       ,2)+q12*fxx(m)       +q13*fxy(m) &
                     +q14*fxz(m)       +q15*fex(m)
                q1f1p=q11*v(n+ninc  ,2)+q12*fxx(m+ninc)  +q13*fxy(m+ninc) &
                     +q14*fxz(m+ninc)  +q15*fex(m+ninc)
                q1f2p=q11*v(n+2*ninc,2)+q12*fxx(m+2*ninc)+q13*fxy(m+2*ninc) &
                     +q14*fxz(m+2*ninc)+q15*fex(m+2*ninc)
!
                q2f1m=q21*v(n-ninc  ,2)+q22*fxx(m-ninc)  +q23*fxy(m-ninc) &
                     +q24*fxz(m-ninc)  +q25*fex(m-ninc)
                q2f  =q21*v(n       ,2)+q22*fxx(m)       +q23*fxy(m) &
                     +q24*fxz(m)       +q25*fex(m)
                q2f1p=q21*v(n+ninc  ,2)+q22*fxx(m+ninc)  +q23*fxy(m+ninc) &
                     +q24*fxz(m+ninc)  +q25*fex(m+ninc)
                q2f2p=q21*v(n+2*ninc,2)+q22*fxx(m+2*ninc)+q23*fxy(m+2*ninc) &
                     +q24*fxz(m+2*ninc)+q25*fex(m+2*ninc)
!
                q3f1m=q31*v(n-ninc  ,2)+q32*fxx(m-ninc)  +q33*fxy(m-ninc) &
                     +q34*fxz(m-ninc)  +q35*fex(m-ninc)
                q3f  =q31*v(n       ,2)+q32*fxx(m)       +q33*fxy(m) &
                     +q34*fxz(m)       +q35*fex(m)
                q3f1p=q31*v(n+ninc  ,2)+q32*fxx(m+ninc)  +q33*fxy(m+ninc) &
                     +q34*fxz(m+ninc)  +q35*fex(m+ninc)
                q3f2p=q31*v(n+2*ninc,2)+q32*fxx(m+2*ninc)+q33*fxy(m+2*ninc) &
                     +q34*fxz(m+2*ninc)+q35*fex(m+2*ninc)
!
                q4f1m=q41*v(n-ninc  ,2)+q42*fxx(m-ninc)  +q43*fxy(m-ninc) &
                     +q44*fxz(m-ninc)  +q45*fex(m-ninc)
                q4f  =q41*v(n       ,2)+q42*fxx(m)       +q43*fxy(m) &
                     +q44*fxz(m)       +q45*fex(m)
                q4f1p=q41*v(n+ninc  ,2)+q42*fxx(m+ninc)  +q43*fxy(m+ninc) &
                     +q44*fxz(m+ninc)  +q45*fex(m+ninc)
                q4f2p=q41*v(n+2*ninc,2)+q42*fxx(m+2*ninc)+q43*fxy(m+2*ninc) &
                     +q44*fxz(m+2*ninc)+q45*fex(m+2*ninc)
!
                q5f1m=q51*v(n-ninc  ,2)+q52*fxx(m-ninc)  +q53*fxy(m-ninc) &
                     +q54*fxz(m-ninc)  +q55*fex(m-ninc)
                q5f  =q51*v(n       ,2)+q52*fxx(m)       +q53*fxy(m) &
                     +q54*fxz(m)       +q55*fex(m)
                q5f1p=q51*v(n+ninc  ,2)+q52*fxx(m+ninc)  +q53*fxy(m+ninc) &
                     +q54*fxz(m+ninc)  +q55*fex(m+ninc)
                q5f2p=q51*v(n+2*ninc,2)+q52*fxx(m+2*ninc)+q53*fxy(m+2*ninc) &
                     +q54*fxz(m+2*ninc)+q55*fex(m+2*ninc)
!        coefficients Crj en maillage irregulier
                c00=cmuj2(m1)*0.5
                c01=cmuj1(m1)*0.5
                c10=-cmuj2(m)*0.5
                c11=1.-c10
                c21=-cmuj1(m+2*ninc)*0.5
                c20=1.-c21
!        calcul des flux d'ordre 2 sur les 2 stencils
                f11=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q1f  *c00 +q1f1p*c01)
                f12=0.5*(1.+sign(1.,v1))*(q1f  *c00 +q1f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q1f1p*c11 +q1f2p*c10)
!
                f21=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q2f  *c00 +q2f1p*c01)
                f22=0.5*(1.+sign(1.,v1))*(q2f  *c00 +q2f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q2f1p*c11 +q2f2p*c10)
!
                f31=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q3f  *c00 +q3f1p*c01)
                f32=0.5*(1.+sign(1.,v1))*(q3f  *c00 +q3f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q3f1p*c11 +q3f2p*c10)
!
                f41=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11) &
                     +0.5*(1.-sign(1.,v4))*(q4f  *c00 +q4f1p*c01)
                f42=0.5*(1.+sign(1.,v4))*(q4f  *c00 +q4f1p*c01) &
                     +0.5*(1.-sign(1.,v4))*(q4f1p*c11 +q4f2p*c10)
!
                f51=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11) &
                     +0.5*(1.-sign(1.,v5))*(q5f  *c00 +q5f1p*c01)
                f52=0.5*(1.+sign(1.,v5))*(q5f  *c00 +q5f1p*c01) &
                     +0.5*(1.-sign(1.,v5))*(q5f1p*c11 +q5f2p*c10)
!        calcul des senseurs beta (au carre)
!         d1c1=(q1f-q1f1m)**2    !maillage regulier
!         d1c2=(q1f1p-q1f)**2
!         d1c3=(q1f2p-q1f1p)**2
                d1c1=4.*((q1f-q1f1m)*c10)**2
                d1c2=4.*((q1f1p-q1f)*c01)**2
                d1c3=4.*((q1f2p-q1f1p)*c01*c21/c00)**2
                d2c1=4.*((q2f-q2f1m)*c10)**2
                d2c2=4.*((q2f1p-q2f)*c01)**2
                d2c3=4.*((q2f2p-q2f1p)*c01*c21/c00)**2
                d3c1=4.*((q3f-q3f1m)*c10)**2
                d3c2=4.*((q3f1p-q3f)*c01)**2
                d3c3=4.*((q3f2p-q3f1p)*c01*c21/c00)**2
                d4c1=4.*((q4f-q4f1m)*c10)**2
                d4c2=4.*((q4f1p-q4f)*c01)**2
                d4c3=4.*((q4f2p-q4f1p)*c01*c21/c00)**2
                d5c1=4.*((q5f-q5f1m)*c10)**2
                d5c2=4.*((q5f1p-q5f)*c01)**2
                d5c3=4.*((q5f2p-q5f1p)*c01*c21/c00)**2
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*d1c1 &
                     +0.5*(1.-sign(1.,v1))*d1c2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*d1c2 &
                     +0.5*(1.-sign(1.,v1))*d1c3+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*d2c1 &
                     +0.5*(1.-sign(1.,v1))*d2c2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*d2c2 &
                     +0.5*(1.-sign(1.,v1))*d2c3+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*d3c1 &
                     +0.5*(1.-sign(1.,v1))*d3c2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*d3c2 &
                     +0.5*(1.-sign(1.,v1))*d3c3+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*d4c1 &
                     +0.5*(1.-sign(1.,v4))*d4c2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*d4c2 &
                     +0.5*(1.-sign(1.,v4))*d4c3+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*d5c1 &
                     +0.5*(1.-sign(1.,v5))*d5c2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*d5c2 &
                     +0.5*(1.-sign(1.,v5))*d5c3+eps)**iexp
!        coefficients gamma en maillage irregulier
                g1p=cmuj2(m1)*cvj(m1)/(cmuj1(m)*cvj(m)+cmuj2(m)*cvj(m)+cmuj2(m1)*cvj(m1))
                g2p=1.-g1p
                g2m=cmuj1(m1)*cvj(m1)/(cmuj1(m1)*cvj(m1)+cmuj2(m1)*cvj(m1)+cmuj2(m+2*ninc)*cvj(m+2*ninc))
                g1m=1.-g2m
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
                ww11m=w11*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w11+w11**2)/(g1p**2+w11*(1.-2.*g1p)) &
                     +w11*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w11+w11**2)/(g1m**2+w11*(1.-2.*g1m))
                ww21m=w21*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w21+w21**2)/(g2p**2+w21*(1.-2.*g2p)) &
                     +w21*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w21+w21**2)/(g2m**2+w21*(1.-2.*g2m))
                swm=ww11m+ww21m
                w11=ww11m/swm
                w21=ww21m/swm
!
                ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
                ww12m=w12*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w12+w12**2)/(g1p**2+w12*(1.-2.*g1p)) &
                     +w12*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w12+w12**2)/(g1m**2+w12*(1.-2.*g1m))
                ww22m=w22*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w22+w22**2)/(g2p**2+w22*(1.-2.*g2p)) &
                     +w22*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w22+w22**2)/(g2m**2+w22*(1.-2.*g2m))
                swm=ww12m+ww22m
                w12=ww12m/swm
                w22=ww22m/swm
!
                ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
                ww13m=w13*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w13+w13**2)/(g1p**2+w13*(1.-2.*g1p)) &
                     +w13*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w13+w13**2)/(g1m**2+w13*(1.-2.*g1m))
                ww23m=w23*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w23+w23**2)/(g2p**2+w23*(1.-2.*g2p)) &
                     +w23*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w23+w23**2)/(g2m**2+w23*(1.-2.*g2m))
                swm=ww13m+ww23m
                w13=ww13m/swm
                w23=ww23m/swm
!
                ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                     +0.5*(1.-sign(1.,v4))*(g1m/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                     +0.5*(1.-sign(1.,v4))*(g2m/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
                ww14m=w14*0.5*(1.+sign(1.,v4))*(g1p+g1p**2-3.*g1p*w14+w14**2)/(g1p**2+w14*(1.-2.*g1p)) &
                     +w14*0.5*(1.-sign(1.,v4))*(g1m+g1m**2-3.*g1m*w14+w14**2)/(g1m**2+w14*(1.-2.*g1m))
                ww24m=w24*0.5*(1.+sign(1.,v4))*(g2p+g2p**2-3.*g2p*w24+w24**2)/(g2p**2+w24*(1.-2.*g2p)) &
                     +w24*0.5*(1.-sign(1.,v4))*(g2m+g2m**2-3.*g2m*w24+w24**2)/(g2m**2+w24*(1.-2.*g2m))
                swm=ww14m+ww24m
                w14=ww14m/swm
                w24=ww24m/swm
!
                ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                     +0.5*(1.-sign(1.,v5))*(g1m/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                     +0.5*(1.-sign(1.,v5))*(g2m/beta52)
                sw=ww15+ww25
                w15=ww15/sw
                w25=ww25/sw
                ww15m=w15*0.5*(1.+sign(1.,v5))*(g1p+g1p**2-3.*g1p*w15+w15**2)/(g1p**2+w15*(1.-2.*g1p)) &
                     +w15*0.5*(1.-sign(1.,v5))*(g1m+g1m**2-3.*g1m*w15+w15**2)/(g1m**2+w15*(1.-2.*g1m))
                ww25m=w25*0.5*(1.+sign(1.,v5))*(g2p+g2p**2-3.*g2p*w25+w25**2)/(g2p**2+w25*(1.-2.*g2p)) &
                     +w25*0.5*(1.-sign(1.,v5))*(g2m+g2m**2-3.*g2m*w25+w25**2)/(g2m**2+w25*(1.-2.*g2m))
                swm=ww15m+ww25m
                w15=ww15m/swm
                w25=ww25m/swm
!        calcul des flux convectifs projetes
                fc1=w11*f11+w21*f12
                fc2=w12*f21+w22*f22
                fc3=w13*f31+w23*f32
                fc4=w14*f41+w24*f42
                fc5=w15*f51+w25*f52
!        produit avec matrice P pour retour dans l'espace physique
                f1=fc1*p11+fc2*p12+fc3*p13+fc4*p14+fc5*p15
                f2=fc1*p21+fc2*p22+fc3*p23+fc4*p24+fc5*p25
                f3=fc1*p31+fc2*p32+fc3*p33+fc4*p34+fc5*p35
                f5=fc1*p51+fc2*p52+fc3*p53+fc4*p54+fc5*p55
!-----------------------------------------------------------------
!        produit de Q avec les flux Euler aux points (i-1) a (i+2)
                q1f1m=q11*v(n-ninc  ,3)+q12*fxy(m-ninc)  +q13*fyy(m-ninc) &
                     +q14*fyz(m-ninc)  +q15*fey(m-ninc)
                q1f  =q11*v(n       ,3)+q12*fxy(m)       +q13*fyy(m) &
                     +q14*fyz(m)       +q15*fey(m)
                q1f1p=q11*v(n+ninc  ,3)+q12*fxy(m+ninc)  +q13*fyy(m+ninc) &
                     +q14*fyz(m+ninc)  +q15*fey(m+ninc)
                q1f2p=q11*v(n+2*ninc,3)+q12*fxy(m+2*ninc)+q13*fyy(m+2*ninc) &
                     +q14*fyz(m+2*ninc)+q15*fey(m+2*ninc)
!
                q2f1m=q21*v(n-ninc  ,3)+q22*fxy(m-ninc)  +q23*fyy(m-ninc) &
                     +q24*fyz(m-ninc)  +q25*fey(m-ninc)
                q2f  =q21*v(n       ,3)+q22*fxy(m)       +q23*fyy(m) &
                     +q24*fyz(m)       +q25*fey(m)
                q2f1p=q21*v(n+ninc  ,3)+q22*fxy(m+ninc)  +q23*fyy(m+ninc) &
                     +q24*fyz(m+ninc)  +q25*fey(m+ninc)
                q2f2p=q21*v(n+2*ninc,3)+q22*fxy(m+2*ninc)+q23*fyy(m+2*ninc) &
                     +q24*fyz(m+2*ninc)+q25*fey(m+2*ninc)
!
                q3f1m=q31*v(n-ninc  ,3)+q32*fxy(m-ninc)  +q33*fyy(m-ninc) &
                     +q34*fyz(m-ninc)  +q35*fey(m-ninc)
                q3f  =q31*v(n       ,3)+q32*fxy(m)       +q33*fyy(m) &
                     +q34*fyz(m)       +q35*fey(m)
                q3f1p=q31*v(n+ninc  ,3)+q32*fxy(m+ninc)  +q33*fyy(m+ninc) &
                     +q34*fyz(m+ninc)  +q35*fey(m+ninc)
                q3f2p=q31*v(n+2*ninc,3)+q32*fxy(m+2*ninc)+q33*fyy(m+2*ninc) &
                     +q34*fyz(m+2*ninc)+q35*fey(m+2*ninc)
!
                q4f1m=q41*v(n-ninc  ,3)+q42*fxy(m-ninc)  +q43*fyy(m-ninc) &
                     +q44*fyz(m-ninc)  +q45*fey(m-ninc)
                q4f  =q41*v(n       ,3)+q42*fxy(m)       +q43*fyy(m) &
                     +q44*fyz(m)       +q45*fey(m)
                q4f1p=q41*v(n+ninc  ,3)+q42*fxy(m+ninc)  +q43*fyy(m+ninc) &
                     +q44*fyz(m+ninc)  +q45*fey(m+ninc)
                q4f2p=q41*v(n+2*ninc,3)+q42*fxy(m+2*ninc)+q43*fyy(m+2*ninc) &
                     +q44*fyz(m+2*ninc)+q45*fey(m+2*ninc)
!
                q5f1m=q51*v(n-ninc  ,3)+q52*fxy(m-ninc)  +q53*fyy(m-ninc) &
                     +q54*fyz(m-ninc)  +q55*fey(m-ninc)
                q5f  =q51*v(n       ,3)+q52*fxy(m)       +q53*fyy(m) &
                     +q54*fyz(m)       +q55*fey(m)
                q5f1p=q51*v(n+ninc  ,3)+q52*fxy(m+ninc)  +q53*fyy(m+ninc) &
                     +q54*fyz(m+ninc)  +q55*fey(m+ninc)
                q5f2p=q51*v(n+2*ninc,3)+q52*fxy(m+2*ninc)+q53*fyy(m+2*ninc) &
                     +q54*fyz(m+2*ninc)+q55*fey(m+2*ninc)
!        calcul des flux d'ordre 2 sur les 2 stencils
                g11=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q1f  *c00 +q1f1p*c01)
                g12=0.5*(1.+sign(1.,v1))*(q1f  *c00 +q1f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q1f1p*c11 +q1f2p*c10)
!
                g21=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q2f  *c00 +q2f1p*c01)
                g22=0.5*(1.+sign(1.,v1))*(q2f  *c00 +q2f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q2f1p*c11 +q2f2p*c10)
!
                g31=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11) &
                     +0.5*(1.-sign(1.,v1))*(q3f  *c00 +q3f1p*c01)
                g32=0.5*(1.+sign(1.,v1))*(q3f  *c00 +q3f1p*c01) &
                     +0.5*(1.-sign(1.,v1))*(q3f1p*c11 +q3f2p*c10)
!
                g41=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11) &
                     +0.5*(1.-sign(1.,v4))*(q4f  *c00 +q4f1p*c01)
                g42=0.5*(1.+sign(1.,v4))*(q4f  *c00 +q4f1p*c01) &
                     +0.5*(1.-sign(1.,v4))*(q4f1p*c11 +q4f2p*c10)
!
                g51=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11) &
                     +0.5*(1.-sign(1.,v5))*(q5f  *c00 +q5f1p*c01)
                g52=0.5*(1.+sign(1.,v5))*(q5f  *c00 +q5f1p*c01) &
                     +0.5*(1.-sign(1.,v5))*(q5f1p*c11 +q5f2p*c10)
!        calcul des senseurs beta (au carre)
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*d1c1 &
                     +0.5*(1.-sign(1.,v1))*d1c2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*d1c2 &
                     +0.5*(1.-sign(1.,v1))*d1c3+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*d2c1 &
                     +0.5*(1.-sign(1.,v1))*d2c2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*d2c2 &
                     +0.5*(1.-sign(1.,v1))*d2c3+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*d3c1 &
                     +0.5*(1.-sign(1.,v1))*d3c2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*d3c2 &
                     +0.5*(1.-sign(1.,v1))*d3c3+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*d4c1 &
                     +0.5*(1.-sign(1.,v4))*d4c2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*d4c2 &
                     +0.5*(1.-sign(1.,v4))*d4c3+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*d5c1 &
                     +0.5*(1.-sign(1.,v5))*d5c2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*d5c2 &
                     +0.5*(1.-sign(1.,v5))*d5c3+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
                ww11m=w11*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w11+w11**2)/(g1p**2+w11*(1.-2.*g1p)) &
                     +w11*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w11+w11**2)/(g1m**2+w11*(1.-2.*g1m))
                ww21m=w21*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w21+w21**2)/(g2p**2+w21*(1.-2.*g2p)) &
                     +w21*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w21+w21**2)/(g2m**2+w21*(1.-2.*g2m))
                swm=ww11m+ww21m
                w11=ww11m/swm
                w21=ww21m/swm
!
                ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
                ww12m=w12*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w12+w12**2)/(g1p**2+w12*(1.-2.*g1p)) &
                     +w12*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w12+w12**2)/(g1m**2+w12*(1.-2.*g1m))
                ww22m=w22*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w22+w22**2)/(g2p**2+w22*(1.-2.*g2p)) &
                     +w22*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w22+w22**2)/(g2m**2+w22*(1.-2.*g2m))
                swm=ww12m+ww22m
                w12=ww12m/swm
                w22=ww22m/swm
!
                ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                     +0.5*(1.-sign(1.,v1))*(g1m/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                     +0.5*(1.-sign(1.,v1))*(g2m/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
                ww13m=w13*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w13+w13**2)/(g1p**2+w13*(1.-2.*g1p)) &
                     +w13*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w13+w13**2)/(g1m**2+w13*(1.-2.*g1m))
                ww23m=w23*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w23+w23**2)/(g2p**2+w23*(1.-2.*g2p)) &
                     +w23*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w23+w23**2)/(g2m**2+w23*(1.-2.*g2m))
                swm=ww13m+ww23m
                w13=ww13m/swm
                w23=ww23m/swm
!
                ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                     +0.5*(1.-sign(1.,v4))*(g1m/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                     +0.5*(1.-sign(1.,v4))*(g2m/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
                ww14m=w14*0.5*(1.+sign(1.,v4))*(g1p+g1p**2-3.*g1p*w14+w14**2)/(g1p**2+w14*(1.-2.*g1p)) &
                     +w14*0.5*(1.-sign(1.,v4))*(g1m+g1m**2-3.*g1m*w14+w14**2)/(g1m**2+w14*(1.-2.*g1m))
                ww24m=w24*0.5*(1.+sign(1.,v4))*(g2p+g2p**2-3.*g2p*w24+w24**2)/(g2p**2+w24*(1.-2.*g2p)) &
                     +w24*0.5*(1.-sign(1.,v4))*(g2m+g2m**2-3.*g2m*w24+w24**2)/(g2m**2+w24*(1.-2.*g2m))
                swm=ww14m+ww24m
                w14=ww14m/swm
                w24=ww24m/swm
!
                ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                     +0.5*(1.-sign(1.,v5))*(g1m/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                     +0.5*(1.-sign(1.,v5))*(g2m/beta52)
                sw=ww15+ww25
                w15=ww15/sw
                w25=ww25/sw
                ww15m=w15*0.5*(1.+sign(1.,v5))*(g1p+g1p**2-3.*g1p*w15+w15**2)/(g1p**2+w15*(1.-2.*g1p)) &
                     +w15*0.5*(1.-sign(1.,v5))*(g1m+g1m**2-3.*g1m*w15+w15**2)/(g1m**2+w15*(1.-2.*g1m))
                ww25m=w25*0.5*(1.+sign(1.,v5))*(g2p+g2p**2-3.*g2p*w25+w25**2)/(g2p**2+w25*(1.-2.*g2p)) &
                     +w25*0.5*(1.-sign(1.,v5))*(g2m+g2m**2-3.*g2m*w25+w25**2)/(g2m**2+w25*(1.-2.*g2m))
                swm=ww15m+ww25m
                w15=ww15m/swm
                w25=ww25m/swm
!        calcul des flux convectifs projetes
                gc1=w11*g11+w21*g12
                gc2=w12*g21+w22*g22
                gc3=w13*g31+w23*g32
                gc4=w14*g41+w24*g42
                gc5=w15*g51+w25*g52
!        produit avec matrice P pour retour dans l'espace physique
                g1=gc1*p11+gc2*p12+gc3*p13+gc4*p14+gc5*p15
                g2=gc1*p21+gc2*p22+gc3*p23+gc4*p24+gc5*p25
                g3=gc1*p31+gc2*p32+gc3*p33+gc4*p34+gc5*p35
                g5=gc1*p51+gc2*p52+gc3*p53+gc4*p54+gc5*p55
!        calcul du flux numerique et bilan de flux
                dg1=f1*sn(m1,kdir,1)+g1*sn(m1,kdir,2)
                dg2=f2*sn(m1,kdir,1)+g2*sn(m1,kdir,2)
                dg3=f3*sn(m1,kdir,1)+g3*sn(m1,kdir,2)
                dg5=f5*sn(m1,kdir,1)+g5*sn(m1,kdir,2)
!        correction point sonique
!         vnl=ul*nx+vl*ny+wl*nz
!         vnr=ur*nx+vr*ny+wr*nz
!         if((vnl-al)*(vnl-ar).lt.0.) then
!          rhoel=v(n,1 )*hl-(ps(n) -pinfl)
!          rhoer=v(n1,1)*hr-(ps(n1)-pinfl)
!          rm=cnds*max(abs(vnl+al),abs(vnr+ar))
!          dj5=rm*(rhoer-rhoel)
!!          dw1=v(n1,1)   -v(n,1)
!!          dw2=v(n1,1)*ur-v(n,1)*ul
!!          dw3=v(n1,1)*vr-v(n,1)*vl
!!          dw4=v(n1,1)*wr-v(n,1)*wl
!!          dw5=rhoer     -rhoel
!!          dj5=(p51*v1*q11+p52*v1*q21+p53*v1*q31+p54*v4*q41+p55*v5*q51)*dw1 &
!!             +(p51*v1*q12+p52*v1*q22+p53*v1*q32+p54*v4*q42+p55*v5*q52)*dw2 &
!!             +(p51*v1*q13+p52*v1*q23+p53*v1*q33+p54*v4*q43+p55*v5*q53)*dw3 &
!!             +(p51*v1*q14+p52*v1*q24+p53*v1*q34+p54*v4*q44+p55*v5*q54)*dw4 &
!!             +(p51*v1*q15+p52*v1*q25+p53*v1*q35+p54*v4*q45+p55*v5*q55)*dw5
!          dfex=cmuj1(m)*(rhoer+ps(n1)-pinfl)*ur  &
!              +cmuj2(m)*(rhoel+ps(n )-pinfl)*ul
!          dfey=cmuj1(m)*(rhoer+ps(n1)-pinfl)*vr  &
!              +cmuj2(m)*(rhoel+ps(n )-pinfl)*vl
!          dg5=0.5*(dfex*sn(m,kdir,1)+dfey*sn(m,kdir,2)-dj5)
!         endif
!        calcul des flux visqueux (multiplies par -2)
                gv2=(cmuj2(m1)*toxx(n)+cmuj1(m1)*toxx(n1))*sn(m1,kdir,1) &
                     +(cmuj2(m1)*toxy(n)+cmuj1(m1)*toxy(n1))*sn(m1,kdir,2)
                gv3=(cmuj2(m1)*toxy(n)+cmuj1(m1)*toxy(n1))*sn(m1,kdir,1) &
                     +(cmuj2(m1)*toyy(n)+cmuj1(m1)*toyy(n1))*sn(m1,kdir,2)
                gv5=(cmuj2(m1)*(toxx(n )*ul+toxy(n )*vl+qcx(n )) &
                     +cmuj1(m1)*(toxx(n1)*ur+toxy(n1)*vr+qcx(n1)))*sn(m1,kdir,1) &
                     +(cmuj2(m1)*(toxy(n )*ul+toyy(n )*vl+qcy(n )) &
                     +cmuj1(m1)*(toxy(n1)*ur+toyy(n1)*vr+qcy(n1)))*sn(m1,kdir,2)
                u(n1,1)=u(n1,1)-dg1
                u(n1,2)=u(n1,2)-dg2+0.5*gv2
                u(n1,3)=u(n1,3)-dg3+0.5*gv3
                u(n1,5)=u(n1,5)-dg5+0.5*gv5
                u(n,1)=u(n,1)+dg1
                u(n,2)=u(n,2)+dg2-0.5*gv2
                u(n,3)=u(n,3)+dg3-0.5*gv3
                u(n,5)=u(n,5)+dg5-0.5*gv5
             enddo
          enddo
       enddo

    endif
!
!-----traitement des bords------------------------------------------
!
    kdir=1
    ninc=nci
!
    do k=k1,k2m1
       ind1 = indc(i1,j1  ,k)
       ind2 = indc(i1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
!       flux a la facette frontiere
          f1=v(n-ninc,2)*sn(m,kdir,1) &
               + v(n-ninc,3)*sn(m,kdir,2)
          f2=(fxx(m-ninc)-toxx(n-ninc))*sn(m,kdir,1) &
               +(fxy(m-ninc)-toxy(n-ninc))*sn(m,kdir,2)
          f3=(fxy(m-ninc)-toxy(n-ninc))*sn(m,kdir,1) &
               +(fyy(m-ninc)-toyy(n-ninc))*sn(m,kdir,2)
          f5=(fex(m-ninc)-(toxx(n-ninc)*v(n-ninc,2)+toxy(n-ninc)* &
               v(n-ninc,3)+ toxz(n-ninc)*v(n-ninc,4))/v(n-ninc,1) &
               -qcx(n-ninc))*sn(m,kdir,1) &
               +(fey(m-ninc)-(toxy(n-ninc)*v(n-ninc,2)+toyy(n-ninc)* &
               v(n-ninc,3)+ toyz(n-ninc)*v(n-ninc,4))/v(n-ninc,1) &
               -qcy(n-ninc))*sn(m,kdir,2)
          u(n,1)=u(n,1)-f1
          u(n,2)=u(n,2)-f2
          u(n,3)=u(n,3)-f3
          u(n,5)=u(n,5)-f5
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i2m1,j1  ,k)
       ind2 = indc(i2m1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          m1=m+ninc
          n1=n+ninc
!       flux a la facette frontiere
          f1=v(n1,2)*sn(m1,kdir,1) &
               + v(n1,3)*sn(m1,kdir,2)
          f2=(fxx(m1)-toxx(n1))*sn(m1,kdir,1) &
               +(fxy(m1)-toxy(n1))*sn(m1,kdir,2)
          f3=(fxy(m1)-toxy(n1))*sn(m1,kdir,1) &
               +(fyy(m1)-toyy(n1))*sn(m1,kdir,2)
          f5=(fex(m1)-(toxx(n1)*v(n1,2)+toxy(n1)*v(n1,3) &
               + toxz(n1)*v(n1,4))/v(n1,1)-qcx(n1))*sn(m1,kdir,1) &
               +(fey(m1)-(toxy(n1)*v(n1,2)+toyy(n1)*v(n1,3) &
               + toyz(n1)*v(n1,4))/v(n1,1)-qcy(n1))*sn(m1,kdir,2)
          u(n,1)=u(n,1)+f1
          u(n,2)=u(n,2)+f2
          u(n,3)=u(n,3)+f3
          u(n,5)=u(n,5)+f5
       enddo
    enddo
!
    kdir=2
    ninc=ncj
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j1,k)
       ind2 = indc(i2m1,j1,k)
!$OMP SIMD
       do n=ind1,ind2
          m=n-n0c
!       flux a la facette frontiere
          g1=v(n-ninc,2)*sn(m,kdir,1) &
               + v(n-ninc,3)*sn(m,kdir,2)
          g2=(fxx(m-ninc)-toxx(n-ninc))*sn(m,kdir,1) &
               +(fxy(m-ninc)-toxy(n-ninc))*sn(m,kdir,2)
          g3=(fxy(m-ninc)-toxy(n-ninc))*sn(m,kdir,1) &
               +(fyy(m-ninc)-toyy(n-ninc))*sn(m,kdir,2)
          g5=(fex(m-ninc)-(toxx(n-ninc)*v(n-ninc,2)+toxy(n-ninc)* &
               v(n-ninc,3)+ toxz(n-ninc)*v(n-ninc,4))/v(n-ninc,1) &
               -qcx(n-ninc))*sn(m,kdir,1) &
               +(fey(m-ninc)-(toxy(n-ninc)*v(n-ninc,2)+toyy(n-ninc)* &
               v(n-ninc,3)+ toyz(n-ninc)*v(n-ninc,4))/v(n-ninc,1) &
               -qcy(n-ninc))*sn(m,kdir,2)
          u(n,1)=u(n,1)-g1
          u(n,2)=u(n,2)-g2
          u(n,3)=u(n,3)-g3
          u(n,5)=u(n,5)-g5
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j2m1,k)
       ind2 = indc(i2m1,j2m1,k)
!$OMP SIMD
       do n=ind1,ind2
          m=n-n0c
          m1=m+ninc
          n1=n+ninc
!       flux a la facette frontiere
          g1=v(n1,2)*sn(m1,kdir,1) &
               + v(n1,3)*sn(m1,kdir,2)
          g2=(fxx(m1)-toxx(n1))*sn(m1,kdir,1) &
               +(fxy(m1)-toxy(n1))*sn(m1,kdir,2)
          g3=(fxy(m1)-toxy(n1))*sn(m1,kdir,1) &
               +(fyy(m1)-toyy(n1))*sn(m1,kdir,2)
          g5=(fex(m1)-(toxx(n1)*v(n1,2)+toxy(n1)*v(n1,3) &
               + toxz(n1)*v(n1,4))/v(n1,1)-qcx(n1))*sn(m1,kdir,1) &
               +(fey(m1)-(toxy(n1)*v(n1,2)+toyy(n1)*v(n1,3) &
               + toyz(n1)*v(n1,4))/v(n1,1)-qcy(n1))*sn(m1,kdir,2)
          u(n,1)=u(n,1)+g1
          u(n,2)=u(n,2)+g2
          u(n,3)=u(n,3)+g3
          u(n,5)=u(n,5)+g5
       enddo
    enddo
!
    if(isortie.eq.1) then
       write(6,'("===>sch_weno3: flux direction i")')
       k=1
       i=13
       do j=j1,j2
          n=indc(i,j,k)
          m=n-n0c
       enddo
!
       write(6,'("===>sch_weno3: flux direction j")')
       k=1
       i=13
       do j=j1,j2
          n=indc(i,j,k)
          m=n-n0c
       enddo
!
       write(6,'("===>sch_weno3: increment explicite")')
       k=1
!       i=13
       i=160
       do j=j1,j2m1
          n=indc(i,j,k)
          m=n-n0c
          write(6,'(i4,i6,4(1pe12.4))') &
               j,n,u(n,1),u(n,2),u(n,4),u(n,5)
       enddo
    endif
!
!-----calcul de la forcing function pour le multigrille--------------
!
    if(ityprk.ne.0) then
       do k=k1,k2m1
          do j=j1,j2m1
             ind1=indc(i1,j,k)
             ind2=indc(i2m1,j,k)
!$OMP SIMD
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
  end subroutine sch_weno3pond
end module mod_sch_weno3pond
