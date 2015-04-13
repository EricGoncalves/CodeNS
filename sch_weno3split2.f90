module mod_sch_weno3split2
  implicit none
contains
  subroutine sch_weno3split2( &
       lm,ityprk, &
       u,v,ff, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       equat, &
       sn,lgsnlt, &
       fxx,fyy,fzz,fxy,fxz,fyz,fex,fey,fez, &
       ps)
!
!******************************************************************
!
!_DA  DATE_C : decembre 2004 - Eric Goncalves / LEGI
!
!     ACT
!_A   Schema WENO de Jiang&Chu, ordre 3 en maillage regulier.
!_A   Formulation 2D avec schema de Roe.
!_A   Flux splitting de Lax-Friedrichs du flux de l'energie pour probleme sonique.
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
    double precision ::               beta41,              beta42,              beta51,              beta52,             betam51
    double precision ::              betam52,             betap51,             betap52,                 c00,                 c01
    double precision ::                  c10,                 c11,                cnds,                 df1,                 df2
    double precision ::                  df3,                 df4,                 df5,                 dg1,                 dg2
    double precision ::                  dg3,                 dg4,                 dg5,                 eps,                  f1
    double precision ::                  f11,                 f12,                  f2,                 f21,                 f22
    double precision ::                   f3,                 f31,                 f32,                  f4,                 f41
    double precision ::                  f42,                  f5,                 f51,                 f52,                 fc1
    double precision ::                  fc2,                 fc3,                 fc4,                 fc5,                fcm5
    double precision ::                 fcp5,           fex(ip00),           fey(ip00),           fez(ip00),       ff(ip11,ip60)
    double precision ::                 fm51,                fm52,                fp51,                fp52,                 fv2
    double precision ::                  fv3,                 fv4,                 fv5,           fxx(ip00),           fxy(ip00)
    double precision ::            fxz(ip00),           fyy(ip00),           fyz(ip00),           fzz(ip00),                  g1
    double precision ::                  g11,                 g12,                  g2,                 g21,                 g22
    double precision ::                   g3,                 g31,                 g32,                  g4,                 g41
    double precision ::                  g42,                  g5,                 g51,                 g52,                 ga1
    double precision ::                  ga2,                 gc1,                 gc2,                 gc3,                 gc4
    double precision ::                  gc5,                  gd,                 gd1,                 gd2,                gm51
    double precision ::                 gm52,                gmc5,                gp51,                gp52,                gpc5
    double precision ::                  gv2,                 gv3,                 gv4,                 gv5,                  hl
    double precision ::                   hm,                  hr,                  nx,                  ny,                  nz
    double precision ::                  p11,                 p12,                 p13,                 p14,                 p15
    double precision ::                  p21,                 p22,                 p23,                 p24,                 p25
    double precision ::                  p31,                 p32,                 p33,                 p34,                 p35
    double precision ::                  p41,                 p42,                 p43,                 p44,                 p45
    double precision ::                  p51,                 p52,                 p53,                 p54,                 p55
    double precision ::             ps(ip11),                 q11,                 q12,                 q13,                 q14
    double precision ::                  q15,                 q1f,               q1f1m,               q1f1p,               q1f2p
    double precision ::                  q21,                 q22,                 q23,                 q24,                 q25
    double precision ::                  q2f,               q2f1m,               q2f1p,               q2f2p,                 q31
    double precision ::                  q32,                 q33,                 q34,                 q35,                 q3f
    double precision ::                q3f1m,               q3f1p,               q3f2p,                 q41,                 q42
    double precision ::                  q43,                 q44,                 q45,                 q4f,               q4f1m
    double precision ::                q4f1p,               q4f2p,                 q51,                 q52,                 q53
    double precision ::                  q54,                 q55,                 q5f,               q5f1m,               q5f1p
    double precision ::                q5f2p,           qcx(ip12),           qcy(ip12),           qcz(ip12),                qm5f
    double precision ::               qm5f1m,              qm5f1p,              qm5f2p,                qp5f,              qp5f1m
    double precision ::               qp5f1p,              qp5f2p,              rhoami,              rhoiam,                rhom
    double precision ::                rhomi,sn(lgsnlt,nind,ndir),                  sw,                 swm,                swmm
    double precision ::                  swp,                swpm,          toxx(ip12),          toxy(ip12),          toxz(ip12)
    double precision ::           toyy(ip12),          toyz(ip12),          tozz(ip12),        u(ip11,ip60),                  ul
    double precision ::                   um,                  ur,        v(ip11,ip60),                  v1,                  v4
    double precision ::                   v5,               vitm2,                  vl,                  vm,                  vn
    double precision ::                   vr,                 w11,                 w12,                 w13,                 w14
    double precision ::                  w15,                 w21,                 w22,                 w23,                 w24
    double precision ::                  w25,                  wl,                  wm,                wm15,                wm25
    double precision ::                 wp15,                wp25,                  wr,                ww11,               ww11m
    double precision ::                 ww12,               ww12m,                ww13,               ww13m,                ww14
    double precision ::                ww14m,                ww15,               ww15m,                ww21,               ww21m
    double precision ::                 ww22,               ww22m,                ww23,               ww23m,                ww24
    double precision ::                ww24m,                ww25,               ww25m,               wwm15,              wwm15m
    double precision ::                wwm25,              wwm25m,               wwp15,              wwp15m,               wwp25
    double precision ::               wwp25m
!
!-------------------------------------------------------------------
!
    character(len=7 ) :: equat
!$OMP MASTER
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
!     coefficient Crj
    c00=0.5
    c01=0.5
    c10=-0.5
    c11=1.5
!     coefficient gamma
    ga1=1./3.
    ga2=2./3.
!     epsilon petit
    eps=1.e-40
!      eps=1.e-6
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
!        splitting du flux de l'energie
                qp5f1m=q51*0.5*(v(n-ninc,2)+abs(v(n-ninc,2))) &
                     +q52*0.5*(fxx(m-ninc)+v(n-ninc,2)*abs(v(n-ninc,2)/v(n-ninc,1))) &
                     +q53*0.5*(fxy(m-ninc)+v(n-ninc,3)*abs(v(n-ninc,2)/v(n-ninc,1))) &
                     +q54*0.5*(fxz(m-ninc)+v(n-ninc,4)*abs(v(n-ninc,2)/v(n-ninc,1) &
                     +sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1))))) &
                     +q55*0.5*(fex(m-ninc)+v(n-ninc,5)*abs(v(n-ninc,2)/v(n-ninc,1) &
                     -sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1)))))
                qm5f1m=q51*0.5*(v(n-ninc,2)-abs(v(n-ninc,2))) &
                     +q52*0.5*(fxx(m-ninc)-v(n-ninc,2)*abs(v(n-ninc,2)/v(n-ninc,1))) &
                     +q53*0.5*(fxy(m-ninc)-v(n-ninc,3)*abs(v(n-ninc,2)/v(n-ninc,1))) &
                     +q54*0.5*(fxz(m-ninc)-v(n-ninc,4)*abs(v(n-ninc,2)/v(n-ninc,1) &
                     +sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1))))) &
                     +q55*0.5*(fex(m-ninc)-v(n-ninc,5)*abs(v(n-ninc,2)/v(n-ninc,1) &
                     -sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1)))))
!
                qp5f=q51*0.5*(v(n,2)+abs(v(n,2))) &
                     +q52*0.5*(fxx(m)+v(n,2)*abs(v(n,2)/v(n,1))) &
                     +q53*0.5*(fxy(m)+v(n,3)*abs(v(n,2)/v(n,1))) &
                     +q54*0.5*(fxz(m)+v(n,4)*abs(v(n,2)/v(n,1)+sqrt(abs(gam*ps(n)/v(n,1))))) &
                     +q55*0.5*(fex(m)+v(n,5)*abs(v(n,2)/v(n,1)-sqrt(abs(gam*ps(n)/v(n,1)))))
                qm5f=q51*0.5*(v(n,2)-abs(v(n,2))) &
                     +q52*0.5*(fxx(m)-v(n,2)*abs(v(n,2)/v(n,1))) &
                     +q53*0.5*(fxy(m)-v(n,3)*abs(v(n,2)/v(n,1))) &
                     +q54*0.5*(fxz(m)-v(n,4)*abs(v(n,2)/v(n,1)+sqrt(abs(gam*ps(n)/v(n,1))))) &
                     +q55*0.5*(fex(m)-v(n,5)*abs(v(n,2)/v(n,1)-sqrt(abs(gam*ps(n)/v(n,1)))))
!
                qp5f1p=q51*0.5*(v(n+ninc,2)+abs(v(n+ninc,2))) &
                     +q52*0.5*(fxx(m+ninc)+v(n+ninc,2)*abs(v(n+ninc,2)/v(n+ninc,1))) &
                     +q53*0.5*(fxy(m+ninc)+v(n+ninc,3)*abs(v(n+ninc,2)/v(n+ninc,1))) &
                     +q54*0.5*(fxz(m+ninc)+v(n+ninc,4)*abs(v(n+ninc,2)/v(n+ninc,1) &
                     +sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1))))) &
                     +q55*0.5*(fex(m+ninc)+v(n+ninc,5)*abs(v(n+ninc,2)/v(n+ninc,1) &
                     -sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1)))))
                qm5f1p=q51*0.5*(v(n+ninc,2)-abs(v(n+ninc,2))) &
                     +q52*0.5*(fxx(m+ninc)-v(n+ninc,2)*abs(v(n+ninc,2)/v(n+ninc,1))) &
                     +q53*0.5*(fxy(m+ninc)-v(n+ninc,3)*abs(v(n+ninc,2)/v(n+ninc,1))) &
                     +q54*0.5*(fxz(m+ninc)-v(n+ninc,4)*abs(v(n+ninc,2)/v(n+ninc,1) &
                     +sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1))))) &
                     +q55*0.5*(fex(m+ninc)-v(n+ninc,5)*abs(v(n+ninc,2)/v(n+ninc,1) &
                     -sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1)))))
!
                qp5f2p=q51*0.5*(v(n+2*ninc,2)+abs(v(n+2*ninc,2))) &
                     +q52*0.5*(fxx(m+2*ninc)+v(n+2*ninc,2)*abs(v(n+2*ninc,2)/v(n+2*ninc,1))) &
                     +q53*0.5*(fxy(m+2*ninc)+v(n+2*ninc,3)*abs(v(n+2*ninc,2)/v(n+2*ninc,1))) &
                     +q54*0.5*(fxz(m+2*ninc)+v(n+2*ninc,4)*abs(v(n+2*ninc,2)/v(n+2*ninc,1) &
                     +sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1))))) &
                     +q55*0.5*(fex(m+2*ninc)+v(n+2*ninc,5)*abs(v(n+2*ninc,2)/v(n+2*ninc,1) &
                     -sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1)))))
                qm5f2p=q51*0.5*(v(n+2*ninc,2)-abs(v(n+2*ninc,2))) &
                     +q52*0.5*(fxx(m+2*ninc)-v(n+2*ninc,2)*abs(v(n+2*ninc,2)/v(n+2*ninc,1))) &
                     +q53*0.5*(fxy(m+2*ninc)-v(n+2*ninc,3)*abs(v(n+2*ninc,2)/v(n+2*ninc,1))) &
                     +q54*0.5*(fxz(m+2*ninc)-v(n+2*ninc,4)*abs(v(n+2*ninc,2)/v(n+2*ninc,1) &
                     +sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1))))) &
                     +q55*0.5*(fex(m+2*ninc)-v(n+2*ninc,5)*abs(v(n+2*ninc,2)/v(n+2*ninc,1) &
                     -sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1)))))
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
                fp51=qp5f1m*c10 +qp5f  *c11
                fm51=qm5f  *c00 +qm5f1p*c01
                fp52=qp5f  *c00 +qp5f1p*c01
                fm52=qm5f1p*c11 +qm5f2p*c10
!        calcul des senseurs beta (au carre)
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*(q1f-q1f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f1p-q1f)**2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*(q1f1p-q1f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f2p-q1f1p)**2+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*(q2f-q2f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f1p-q2f)**2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*(q2f1p-q2f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f2p-q2f1p)**2+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*(q3f-q3f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f1p-q3f)**2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*(q3f1p-q3f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f2p-q3f1p)**2+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*(q4f-q4f1m)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f1p-q4f)**2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*(q4f1p-q4f)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f2p-q4f1p)**2+eps)**iexp
!
                betap51=((qp5f  -qp5f1m)**2+eps)**iexp
                betam51=((qm5f1p-qm5f  )**2+eps)**iexp
                betap52=((qp5f1p-qp5f  )**2+eps)**iexp
                betam52=((qm5f2p-qm5f1p)**2+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(ga1/beta11) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(ga2/beta12) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
!
                ww12=0.5*(1.+sign(1.,v1))*(ga1/beta21) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(ga2/beta22) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
!
                ww13=0.5*(1.+sign(1.,v1))*(ga1/beta31) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(ga2/beta32) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
!
                ww14=0.5*(1.+sign(1.,v4))*(ga1/beta41) &
                     +0.5*(1.-sign(1.,v4))*(ga2/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(ga2/beta42) &
                     +0.5*(1.-sign(1.,v4))*(ga1/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
!
                wwp15=ga1/betap51
                wwm15=ga2/betam51
                wwp25=ga2/betap52
                wwm25=ga1/betam52
                swp=wwp15+wwp25
                swm=wwm15+wwm25
                wp15=wwp15/swp
                wp25=wwp25/swp
                wm15=wwm15/swm
                wm25=wwm25/swm
!        calcul des flux convectifs projetes
                fc1=w11*f11+w21*f12
                fc2=w12*f21+w22*f22
                fc3=w13*f31+w23*f32
                fc4=w14*f41+w24*f42
                fcp5=wp15*fp51+wp25*fp52
                fcm5=wm15*fm51+wm25*fm52
                fc5=fcp5+fcm5
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
                beta11=(0.5*(1.+sign(1.,v1))*(q1f-q1f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f1p-q1f)**2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*(q1f1p-q1f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f2p-q1f1p)**2+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*(q2f-q2f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f1p-q2f)**2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*(q2f1p-q2f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f2p-q2f1p)**2+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*(q3f-q3f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f1p-q3f)**2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*(q3f1p-q3f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f2p-q3f1p)**2+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*(q4f-q4f1m)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f1p-q4f)**2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*(q4f1p-q4f)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f2p-q4f1p)**2+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*(q5f-q5f1m)**2 &
                     +0.5*(1.-sign(1.,v5))*(q5f1p-q5f)**2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*(q5f1p-q5f)**2 &
                     +0.5*(1.-sign(1.,v5))*(q5f2p-q5f1p)**2+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(ga1/beta11) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(ga2/beta12) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
!
                ww12=0.5*(1.+sign(1.,v1))*(ga1/beta21) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(ga2/beta22) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
!
                ww13=0.5*(1.+sign(1.,v1))*(ga1/beta31) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(ga2/beta32) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
!
                ww14=0.5*(1.+sign(1.,v4))*(ga1/beta41) &
                     +0.5*(1.-sign(1.,v4))*(ga2/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(ga2/beta42) &
                     +0.5*(1.-sign(1.,v4))*(ga1/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
!
                ww15=0.5*(1.+sign(1.,v5))*(ga1/beta51) &
                     +0.5*(1.-sign(1.,v5))*(ga2/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(ga2/beta52) &
                     +0.5*(1.-sign(1.,v5))*(ga1/beta52)
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
                fv2=(toxx(n)+toxx(n1))*sn(m1,kdir,1) &
                     +(toxy(n)+toxy(n1))*sn(m1,kdir,2)
                fv3=(toxy(n)+toxy(n1))*sn(m1,kdir,1) &
                     +(toyy(n)+toyy(n1))*sn(m1,kdir,2)
                fv5=(toxx(n )*ul+toxy(n )*vl+qcx(n ) &
                     +toxx(n1)*ur+toxy(n1)*vr+qcx(n1))*sn(m1,kdir,1) &
                     +(toxy(n )*ul+toyy(n )*vl+qcy(n ) &
                     +toxy(n1)*ur+toyy(n1)*vr+qcy(n1))*sn(m1,kdir,2)
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
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*(q1f-q1f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f1p-q1f)**2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*(q1f1p-q1f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f2p-q1f1p)**2+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*(q2f-q2f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f1p-q2f)**2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*(q2f1p-q2f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f2p-q2f1p)**2+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*(q3f-q3f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f1p-q3f)**2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*(q3f1p-q3f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f2p-q3f1p)**2+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*(q4f-q4f1m)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f1p-q4f)**2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*(q4f1p-q4f)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f2p-q4f1p)**2+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*(q5f-q5f1m)**2 &
                     +0.5*(1.-sign(1.,v5))*(q5f1p-q5f)**2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*(q5f1p-q5f)**2 &
                     +0.5*(1.-sign(1.,v5))*(q5f2p-q5f1p)**2+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(ga1/beta11) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(ga2/beta12) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
!
                ww12=0.5*(1.+sign(1.,v1))*(ga1/beta21) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(ga2/beta22) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
!
                ww13=0.5*(1.+sign(1.,v1))*(ga1/beta31) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(ga2/beta32) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
!
                ww14=0.5*(1.+sign(1.,v4))*(ga1/beta41) &
                     +0.5*(1.-sign(1.,v4))*(ga2/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(ga2/beta42) &
                     +0.5*(1.-sign(1.,v4))*(ga1/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
!
                ww15=0.5*(1.+sign(1.,v5))*(ga1/beta51) &
                     +0.5*(1.-sign(1.,v5))*(ga2/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(ga2/beta52) &
                     +0.5*(1.-sign(1.,v5))*(ga1/beta52)
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
!        splitting du flux de l'energie
                qp5f1m=q51*0.5*(v(n-ninc,3)+abs(v(n-ninc,3))) &
                     +q52*0.5*(fxy(m-ninc)+v(n-ninc,2)*abs(v(n-ninc,3)/v(n-ninc,1))) &
                     +q53*0.5*(fyy(m-ninc)+v(n-ninc,3)*abs(v(n-ninc,3)/v(n-ninc,1))) &
                     +q54*0.5*(fyz(m-ninc)+v(n-ninc,4)*abs(v(n-ninc,3)/v(n-ninc,1) &
                     +sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1))))) &
                     +q55*0.5*(fey(m-ninc)+v(n-ninc,5)*abs(v(n-ninc,3)/v(n-ninc,1) &
                     -sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1)))))
                qm5f1m=q51*0.5*(v(n-ninc,3)-abs(v(n-ninc,3))) &
                     +q52*0.5*(fxy(m-ninc)-v(n-ninc,2)*abs(v(n-ninc,3)/v(n-ninc,1))) &
                     +q53*0.5*(fyy(m-ninc)-v(n-ninc,3)*abs(v(n-ninc,3)/v(n-ninc,1))) &
                     +q54*0.5*(fyz(m-ninc)-v(n-ninc,4)*abs(v(n-ninc,3)/v(n-ninc,1) &
                     +sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1))))) &
                     +q55*0.5*(fey(m-ninc)-v(n-ninc,5)*abs(v(n-ninc,3)/v(n-ninc,1) &
                     -sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1)))))
!
                qp5f=q51*0.5*(v(n,3)+abs(v(n,3))) &
                     +q52*0.5*(fxy(m)+v(n,2)*abs(v(n,3)/v(n,1))) &
                     +q53*0.5*(fyy(m)+v(n,3)*abs(v(n,3)/v(n,1))) &
                     +q54*0.5*(fyz(m)+v(n,4)*abs(v(n,3)/v(n,1)+sqrt(abs(gam*ps(n)/v(n,1))))) &
                     +q55*0.5*(fey(m)+v(n,5)*abs(v(n,3)/v(n,1)-sqrt(abs(gam*ps(n)/v(n,1)))))
                qm5f=q51*0.5*(v(n,3)-abs(v(n,3))) &
                     +q52*0.5*(fxy(m)-v(n,2)*abs(v(n,3)/v(n,1))) &
                     +q53*0.5*(fyy(m)-v(n,3)*abs(v(n,3)/v(n,1))) &
                     +q54*0.5*(fyz(m)-v(n,4)*abs(v(n,3)/v(n,1)+sqrt(abs(gam*ps(n)/v(n,1))))) &
                     +q55*0.5*(fey(m)-v(n,5)*abs(v(n,3)/v(n,1)-sqrt(abs(gam*ps(n)/v(n,1)))))
!
                qp5f1p=q51*0.5*(v(n+ninc,3)+abs(v(n+ninc,3))) &
                     +q52*0.5*(fxy(m+ninc)+v(n+ninc,2)*abs(v(n+ninc,3)/v(n+ninc,1))) &
                     +q53*0.5*(fyy(m+ninc)+v(n+ninc,3)*abs(v(n+ninc,3)/v(n+ninc,1))) &
                     +q54*0.5*(fyz(m+ninc)+v(n+ninc,4)*abs(v(n+ninc,3)/v(n+ninc,1) &
                     +sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1))))) &
                     +q55*0.5*(fey(m+ninc)+v(n+ninc,5)*abs(v(n+ninc,3)/v(n+ninc,1) &
                     -sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1)))))
                qm5f1p=q51*0.5*(v(n+ninc,3)-abs(v(n+ninc,3))) &
                     +q52*0.5*(fxy(m+ninc)-v(n+ninc,2)*abs(v(n+ninc,3)/v(n+ninc,1))) &
                     +q53*0.5*(fyy(m+ninc)-v(n+ninc,3)*abs(v(n+ninc,3)/v(n+ninc,1))) &
                     +q54*0.5*(fyz(m+ninc)-v(n+ninc,4)*abs(v(n+ninc,3)/v(n+ninc,1) &
                     +sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1))))) &
                     +q55*0.5*(fey(m+ninc)-v(n+ninc,5)*abs(v(n+ninc,3)/v(n+ninc,1) &
                     -sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1)))))
!
                qp5f2p=q51*0.5*(v(n+2*ninc,3)+abs(v(n+2*ninc,3))) &
                     +q52*0.5*(fxy(m+2*ninc)+v(n+2*ninc,2)*abs(v(n+2*ninc,3)/v(n+2*ninc,1))) &
                     +q53*0.5*(fyy(m+2*ninc)+v(n+2*ninc,3)*abs(v(n+2*ninc,3)/v(n+2*ninc,1))) &
                     +q54*0.5*(fyz(m+2*ninc)+v(n+2*ninc,4)*abs(v(n+2*ninc,3)/v(n+2*ninc,1) &
                     +sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1))))) &
                     +q55*0.5*(fey(m+2*ninc)+v(n+2*ninc,5)*abs(v(n+2*ninc,3)/v(n+2*ninc,1) &
                     -sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1)))))
                qm5f2p=q51*0.5*(v(n+2*ninc,3)-abs(v(n+2*ninc,3))) &
                     +q52*0.5*(fxy(m+2*ninc)-v(n+2*ninc,2)*abs(v(n+2*ninc,3)/v(n+2*ninc,1))) &
                     +q53*0.5*(fyy(m+2*ninc)-v(n+2*ninc,3)*abs(v(n+2*ninc,3)/v(n+2*ninc,1))) &
                     +q54*0.5*(fyz(m+2*ninc)-v(n+2*ninc,4)*abs(v(n+2*ninc,3)/v(n+2*ninc,1) &
                     +sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1))))) &
                     +q55*0.5*(fey(m+2*ninc)-v(n+2*ninc,5)*abs(v(n+2*ninc,3)/v(n+2*ninc,1) &
                     -sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1)))))
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
                gp51=qp5f1m*c10 +qp5f  *c11
                gm51=qm5f  *c00 +qm5f1p*c01
                gp52=qp5f  *c00 +qp5f1p*c01
                gm52=qm5f1p*c11 +qm5f2p*c10
!        calcul des senseurs beta (au carre)
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*(q1f-q1f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f1p-q1f)**2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*(q1f1p-q1f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f2p-q1f1p)**2+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*(q2f-q2f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f1p-q2f)**2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*(q2f1p-q2f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f2p-q2f1p)**2+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*(q3f-q3f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f1p-q3f)**2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*(q3f1p-q3f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f2p-q3f1p)**2+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*(q4f-q4f1m)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f1p-q4f)**2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*(q4f1p-q4f)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f2p-q4f1p)**2+eps)**iexp
!
                betap51=((qp5f  -qp5f1m)**2+eps)**iexp
                betam51=((qm5f1p-qm5f  )**2+eps)**iexp
                betap52=((qp5f1p-qp5f  )**2+eps)**iexp
                betam52=((qm5f2p-qm5f1p)**2+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(ga1/beta11) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(ga2/beta12) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
!
                ww12=0.5*(1.+sign(1.,v1))*(ga1/beta21) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(ga2/beta22) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
!
                ww13=0.5*(1.+sign(1.,v1))*(ga1/beta31) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(ga2/beta32) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
!
                ww14=0.5*(1.+sign(1.,v4))*(ga1/beta41) &
                     +0.5*(1.-sign(1.,v4))*(ga2/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(ga2/beta42) &
                     +0.5*(1.-sign(1.,v4))*(ga1/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
!
                wwp15=ga1/betap51
                wwm15=ga2/betam51
                wwp25=ga2/betap52
                wwm25=ga1/betam52
                swp=wwp15+wwp25
                swm=wwm15+wwm25
                wp15=wwp15/swp
                wp25=wwp25/swp
                wm15=wwm15/swm
                wm25=wwm25/swm
!        calcul des flux convectifs projetes
                gc1=w11*g11+w21*g12
                gc2=w12*g21+w22*g22
                gc3=w13*g31+w23*g32
                gc4=w14*g41+w24*g42
                gpc5=wp15*gp51+wp25*gp52
                gmc5=wm15*gm51+wm25*gm52
                gc5=gpc5+gmc5
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
                gv2=(toxx(n)+toxx(n1))*sn(m1,kdir,1) &
                     +(toxy(n)+toxy(n1))*sn(m1,kdir,2)
                gv3=(toxy(n)+toxy(n1))*sn(m1,kdir,1) &
                     +(toyy(n)+toyy(n1))*sn(m1,kdir,2)
                gv5=(toxx(n )*ul+toxy(n )*vl+qcx(n ) &
                     +toxx(n1)*ur+toxy(n1)*vr+qcx(n1))*sn(m1,kdir,1) &
                     +(toxy(n )*ul+toyy(n )*vl+qcy(n ) &
                     +toxy(n1)*ur+toyy(n1)*vr+qcy(n1))*sn(m1,kdir,2)
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
!
    elseif(imap.eq.1) then !mapping de henrik
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
!        splitting du flux de l'energie
                qp5f1m=q51*0.5*(v(n-ninc,2)+abs(v(n-ninc,2))) &
                     +q52*0.5*(fxx(m-ninc)+v(n-ninc,2)*abs(v(n-ninc,2)/v(n-ninc,1))) &
                     +q53*0.5*(fxy(m-ninc)+v(n-ninc,3)*abs(v(n-ninc,2)/v(n-ninc,1))) &
                     +q54*0.5*(fxz(m-ninc)+v(n-ninc,4)*abs(v(n-ninc,2)/v(n-ninc,1) &
                     +sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1))))) &
                     +q55*0.5*(fex(m-ninc)+v(n-ninc,5)*abs(v(n-ninc,2)/v(n-ninc,1) &
                     -sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1)))))
                qm5f1m=q51*0.5*(v(n-ninc,2)-abs(v(n-ninc,2))) &
                     +q52*0.5*(fxx(m-ninc)-v(n-ninc,2)*abs(v(n-ninc,2)/v(n-ninc,1))) &
                     +q53*0.5*(fxy(m-ninc)-v(n-ninc,3)*abs(v(n-ninc,2)/v(n-ninc,1))) &
                     +q54*0.5*(fxz(m-ninc)-v(n-ninc,4)*abs(v(n-ninc,2)/v(n-ninc,1) &
                     +sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1))))) &
                     +q55*0.5*(fex(m-ninc)-v(n-ninc,5)*abs(v(n-ninc,2)/v(n-ninc,1) &
                     -sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1)))))
!
                qp5f=q51*0.5*(v(n,2)+abs(v(n,2))) &
                     +q52*0.5*(fxx(m)+v(n,2)*abs(v(n,2)/v(n,1))) &
                     +q53*0.5*(fxy(m)+v(n,3)*abs(v(n,2)/v(n,1))) &
                     +q54*0.5*(fxz(m)+v(n,4)*abs(v(n,2)/v(n,1)+sqrt(abs(gam*ps(n)/v(n,1))))) &
                     +q55*0.5*(fex(m)+v(n,5)*abs(v(n,2)/v(n,1)-sqrt(abs(gam*ps(n)/v(n,1)))))
                qm5f=q51*0.5*(v(n,2)-abs(v(n,2))) &
                     +q52*0.5*(fxx(m)-v(n,2)*abs(v(n,2)/v(n,1))) &
                     +q53*0.5*(fxy(m)-v(n,3)*abs(v(n,2)/v(n,1))) &
                     +q54*0.5*(fxz(m)-v(n,4)*abs(v(n,2)/v(n,1)+sqrt(abs(gam*ps(n)/v(n,1))))) &
                     +q55*0.5*(fex(m)-v(n,5)*abs(v(n,2)/v(n,1)-sqrt(abs(gam*ps(n)/v(n,1)))))
!
                qp5f1p=q51*0.5*(v(n+ninc,2)+abs(v(n+ninc,2))) &
                     +q52*0.5*(fxx(m+ninc)+v(n+ninc,2)*abs(v(n+ninc,2)/v(n+ninc,1))) &
                     +q53*0.5*(fxy(m+ninc)+v(n+ninc,3)*abs(v(n+ninc,2)/v(n+ninc,1))) &
                     +q54*0.5*(fxz(m+ninc)+v(n+ninc,4)*abs(v(n+ninc,2)/v(n+ninc,1) &
                     +sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1))))) &
                     +q55*0.5*(fex(m+ninc)+v(n+ninc,5)*abs(v(n+ninc,2)/v(n+ninc,1) &
                     -sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1)))))
                qm5f1p=q51*0.5*(v(n+ninc,2)-abs(v(n+ninc,2))) &
                     +q52*0.5*(fxx(m+ninc)-v(n+ninc,2)*abs(v(n+ninc,2)/v(n+ninc,1))) &
                     +q53*0.5*(fxy(m+ninc)-v(n+ninc,3)*abs(v(n+ninc,2)/v(n+ninc,1))) &
                     +q54*0.5*(fxz(m+ninc)-v(n+ninc,4)*abs(v(n+ninc,2)/v(n+ninc,1) &
                     +sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1))))) &
                     +q55*0.5*(fex(m+ninc)-v(n+ninc,5)*abs(v(n+ninc,2)/v(n+ninc,1) &
                     -sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1)))))
!
                qp5f2p=q51*0.5*(v(n+2*ninc,2)+abs(v(n+2*ninc,2))) &
                     +q52*0.5*(fxx(m+2*ninc)+v(n+2*ninc,2)*abs(v(n+2*ninc,2)/v(n+2*ninc,1))) &
                     +q53*0.5*(fxy(m+2*ninc)+v(n+2*ninc,3)*abs(v(n+2*ninc,2)/v(n+2*ninc,1))) &
                     +q54*0.5*(fxz(m+2*ninc)+v(n+2*ninc,4)*abs(v(n+2*ninc,2)/v(n+2*ninc,1) &
                     +sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1))))) &
                     +q55*0.5*(fex(m+2*ninc)+v(n+2*ninc,5)*abs(v(n+2*ninc,2)/v(n+2*ninc,1) &
                     -sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1)))))
                qm5f2p=q51*0.5*(v(n+2*ninc,2)-abs(v(n+2*ninc,2))) &
                     +q52*0.5*(fxx(m+2*ninc)-v(n+2*ninc,2)*abs(v(n+2*ninc,2)/v(n+2*ninc,1))) &
                     +q53*0.5*(fxy(m+2*ninc)-v(n+2*ninc,3)*abs(v(n+2*ninc,2)/v(n+2*ninc,1))) &
                     +q54*0.5*(fxz(m+2*ninc)-v(n+2*ninc,4)*abs(v(n+2*ninc,2)/v(n+2*ninc,1) &
                     +sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1))))) &
                     +q55*0.5*(fex(m+2*ninc)-v(n+2*ninc,5)*abs(v(n+2*ninc,2)/v(n+2*ninc,1) &
                     -sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1)))))
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
                fp51=qp5f1m*c10 +qp5f  *c11
                fm51=qm5f  *c00 +qm5f1p*c01
                fp52=qp5f  *c00 +qp5f1p*c01
                fm52=qm5f1p*c11 +qm5f2p*c10
!        calcul des senseurs beta (au carre)
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*(q1f-q1f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f1p-q1f)**2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*(q1f1p-q1f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f2p-q1f1p)**2+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*(q2f-q2f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f1p-q2f)**2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*(q2f1p-q2f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f2p-q2f1p)**2+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*(q3f-q3f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f1p-q3f)**2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*(q3f1p-q3f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f2p-q3f1p)**2+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*(q4f-q4f1m)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f1p-q4f)**2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*(q4f1p-q4f)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f2p-q4f1p)**2+eps)**iexp
!
                betap51=((qp5f  -qp5f1m)**2+eps)**iexp
                betam51=((qm5f1p-qm5f  )**2+eps)**iexp
                betap52=((qp5f1p-qp5f  )**2+eps)**iexp
                betam52=((qm5f2p-qm5f1p)**2+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(ga1/beta11) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(ga2/beta12) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
                ww11m=w11*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w11+w11**2)/(ga1**2+w11*(1.-2.*ga1)) &
                     +w11*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w11+w11**2)/(ga2**2+w11*(1.-2.*ga2))
                ww21m=w21*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w21+w21**2)/(ga2**2+w21*(1.-2.*ga2)) &
                     +w21*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w21+w21**2)/(ga1**2+w21*(1.-2.*ga1))
                swm=ww11m+ww21m
                w11=ww11m/swm
                w21=ww21m/swm
!
                ww12=0.5*(1.+sign(1.,v1))*(ga1/beta21) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(ga2/beta22) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
                ww12m=w12*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w12+w12**2)/(ga1**2+w12*(1.-2.*ga1)) &
                     +w12*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w12+w12**2)/(ga2**2+w12*(1.-2.*ga2))
                ww22m=w22*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w22+w22**2)/(ga2**2+w22*(1.-2.*ga2)) &
                     +w22*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w22+w22**2)/(ga1**2+w22*(1.-2.*ga1))
                swm=ww12m+ww22m
                w12=ww12m/swm
                w22=ww22m/swm
!
                ww13=0.5*(1.+sign(1.,v1))*(ga1/beta31) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(ga2/beta32) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
                ww13m=w13*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w13+w13**2)/(ga1**2+w13*(1.-2.*ga1)) &
                     +w13*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w13+w13**2)/(ga2**2+w13*(1.-2.*ga2))
                ww23m=w23*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w23+w23**2)/(ga2**2+w23*(1.-2.*ga2)) &
                     +w23*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w23+w23**2)/(ga1**2+w23*(1.-2.*ga1))
                swm=ww13m+ww23m
                w13=ww13m/swm
                w23=ww23m/swm
!
                ww14=0.5*(1.+sign(1.,v4))*(ga1/beta41) &
                     +0.5*(1.-sign(1.,v4))*(ga2/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(ga2/beta42) &
                     +0.5*(1.-sign(1.,v4))*(ga1/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
                ww14m=w14*0.5*(1.+sign(1.,v4))*(ga1+ga1**2-3.*ga1*w14+w14**2)/(ga1**2+w14*(1.-2.*ga1)) &
                     +w14*0.5*(1.-sign(1.,v4))*(ga2+ga2**2-3.*ga2*w14+w14**2)/(ga2**2+w14*(1.-2.*ga2))
                ww24m=w24*0.5*(1.+sign(1.,v4))*(ga2+ga2**2-3.*ga2*w24+w24**2)/(ga2**2+w24*(1.-2.*ga2)) &
                     +w24*0.5*(1.-sign(1.,v4))*(ga1+ga1**2-3.*ga1*w24+w24**2)/(ga1**2+w24*(1.-2.*ga1))
                swm=ww14m+ww24m
                w14=ww14m/swm
                w24=ww24m/swm
!
                wwp15=ga1/betap51
                wwm15=ga2/betam51
                wwp25=ga2/betap52
                wwm25=ga1/betam52
                swp=wwp15+wwp25
                swm=wwm15+wwm25
                wp15=wwp15/swp
                wp25=wwp25/swp
                wm15=wwm15/swm
                wm25=wwm25/swm
                wwp15m=wp15*(ga1+ga1**2-3.*ga1*wp15+wp15**2)/(ga1**2+wp15*(1.-2.*ga1))
                wwm15m=wm15*(ga2+ga2**2-3.*ga2*wm15+wm15**2)/(ga2**2+wm15*(1.-2.*ga2))
                wwp25m=wp25*(ga2+ga2**2-3.*ga2*wp25+wp25**2)/(ga2**2+wp25*(1.-2.*ga2))
                wwm25m=wm25*(ga1+ga1**2-3.*ga1*wm25+wm25**2)/(ga1**2+wm25*(1.-2.*ga1))
                swpm=wwp15m+wwp25m
                swmm=wwm15m+wwm25m
                wp15=wwp15m/swpm
                wp25=wwp25m/swpm
                wm15=wwm15m/swmm
                wm25=wwm25m/swmm
!        calcul des flux convectifs projetes
                fc1=w11*f11+w21*f12
                fc2=w12*f21+w22*f22
                fc3=w13*f31+w23*f32
                fc4=w14*f41+w24*f42
                fcp5=wp15*fp51+wp25*fp52
                fcm5=wm15*fm51+wm25*fm52
                fc5=fcp5+fcm5
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
                beta11=(0.5*(1.+sign(1.,v1))*(q1f-q1f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f1p-q1f)**2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*(q1f1p-q1f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f2p-q1f1p)**2+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*(q2f-q2f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f1p-q2f)**2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*(q2f1p-q2f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f2p-q2f1p)**2+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*(q3f-q3f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f1p-q3f)**2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*(q3f1p-q3f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f2p-q3f1p)**2+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*(q4f-q4f1m)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f1p-q4f)**2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*(q4f1p-q4f)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f2p-q4f1p)**2+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*(q5f-q5f1m)**2 &
                     +0.5*(1.-sign(1.,v5))*(q5f1p-q5f)**2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*(q5f1p-q5f)**2 &
                     +0.5*(1.-sign(1.,v5))*(q5f2p-q5f1p)**2+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(ga1/beta11) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(ga2/beta12) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
                ww11m=w11*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w11+w11**2)/(ga1**2+w11*(1.-2.*ga1)) &
                     +w11*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w11+w11**2)/(ga2**2+w11*(1.-2.*ga2))
                ww21m=w21*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w21+w21**2)/(ga2**2+w21*(1.-2.*ga2)) &
                     +w21*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w21+w21**2)/(ga1**2+w21*(1.-2.*ga1))
                swm=ww11m+ww21m
                w11=ww11m/swm
                w21=ww21m/swm
!
                ww12=0.5*(1.+sign(1.,v1))*(ga1/beta21) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(ga2/beta22) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
                ww12m=w12*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w12+w12**2)/(ga1**2+w12*(1.-2.*ga1)) &
                     +w12*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w12+w12**2)/(ga2**2+w12*(1.-2.*ga2))
                ww22m=w22*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w22+w22**2)/(ga2**2+w22*(1.-2.*ga2)) &
                     +w22*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w22+w22**2)/(ga1**2+w22*(1.-2.*ga1))
                swm=ww12m+ww22m
                w12=ww12m/swm
                w22=ww22m/swm
!
                ww13=0.5*(1.+sign(1.,v1))*(ga1/beta31) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(ga2/beta32) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
                ww13m=w13*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w13+w13**2)/(ga1**2+w13*(1.-2.*ga1)) &
                     +w13*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w13+w13**2)/(ga2**2+w13*(1.-2.*ga2))
                ww23m=w23*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w23+w23**2)/(ga2**2+w23*(1.-2.*ga2)) &
                     +w23*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w23+w23**2)/(ga1**2+w23*(1.-2.*ga1))
                swm=ww13m+ww23m
                w13=ww13m/swm
                w23=ww23m/swm
!
                ww14=0.5*(1.+sign(1.,v4))*(ga1/beta41) &
                     +0.5*(1.-sign(1.,v4))*(ga2/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(ga2/beta42) &
                     +0.5*(1.-sign(1.,v4))*(ga1/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
                ww14m=w14*0.5*(1.+sign(1.,v4))*(ga1+ga1**2-3.*ga1*w14+w14**2)/(ga1**2+w14*(1.-2.*ga1)) &
                     +w14*0.5*(1.-sign(1.,v4))*(ga2+ga2**2-3.*ga2*w14+w14**2)/(ga2**2+w14*(1.-2.*ga2))
                ww24m=w24*0.5*(1.+sign(1.,v4))*(ga2+ga2**2-3.*ga2*w24+w24**2)/(ga2**2+w24*(1.-2.*ga2)) &
                     +w24*0.5*(1.-sign(1.,v4))*(ga1+ga1**2-3.*ga1*w24+w24**2)/(ga1**2+w24*(1.-2.*ga1))
                swm=ww14m+ww24m
                w14=ww14m/swm
                w24=ww24m/swm
!
                ww15=0.5*(1.+sign(1.,v5))*(ga1/beta51) &
                     +0.5*(1.-sign(1.,v5))*(ga2/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(ga2/beta52) &
                     +0.5*(1.-sign(1.,v5))*(ga1/beta52)
                sw=ww15+ww25
                w15=ww15/sw
                w25=ww25/sw
                ww15m=w15*0.5*(1.+sign(1.,v5))*(ga1+ga1**2-3.*ga1*w15+w15**2)/(ga1**2+w15*(1.-2.*ga1)) &
                     +w15*0.5*(1.-sign(1.,v5))*(ga2+ga2**2-3.*ga2*w15+w15**2)/(ga2**2+w15*(1.-2.*ga2))
                ww25m=w25*0.5*(1.+sign(1.,v5))*(ga2+ga2**2-3.*ga2*w25+w25**2)/(ga2**2+w25*(1.-2.*ga2)) &
                     +w25*0.5*(1.-sign(1.,v5))*(ga1+ga1**2-3.*ga1*w25+w25**2)/(ga1**2+w25*(1.-2.*ga1))
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
!        calcul des flux visqueux (multiplies par -2)
                fv2=(toxx(n)+toxx(n1))*sn(m1,kdir,1) &
                     +(toxy(n)+toxy(n1))*sn(m1,kdir,2)
                fv3=(toxy(n)+toxy(n1))*sn(m1,kdir,1) &
                     +(toyy(n)+toyy(n1))*sn(m1,kdir,2)
                fv5=(toxx(n )*ul+toxy(n )*vl+qcx(n ) &
                     +toxx(n1)*ur+toxy(n1)*vr+qcx(n1))*sn(m1,kdir,1) &
                     +(toxy(n )*ul+toyy(n )*vl+qcy(n ) &
                     +toxy(n1)*ur+toyy(n1)*vr+qcy(n1))*sn(m1,kdir,2)
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
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*(q1f-q1f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f1p-q1f)**2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*(q1f1p-q1f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f2p-q1f1p)**2+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*(q2f-q2f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f1p-q2f)**2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*(q2f1p-q2f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f2p-q2f1p)**2+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*(q3f-q3f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f1p-q3f)**2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*(q3f1p-q3f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f2p-q3f1p)**2+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*(q4f-q4f1m)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f1p-q4f)**2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*(q4f1p-q4f)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f2p-q4f1p)**2+eps)**iexp
!
                beta51=(0.5*(1.+sign(1.,v5))*(q5f-q5f1m)**2 &
                     +0.5*(1.-sign(1.,v5))*(q5f1p-q5f)**2+eps)**iexp
                beta52=(0.5*(1.+sign(1.,v5))*(q5f1p-q5f)**2 &
                     +0.5*(1.-sign(1.,v5))*(q5f2p-q5f1p)**2+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(ga1/beta11) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(ga2/beta12) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
                ww11m=w11*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w11+w11**2)/(ga1**2+w11*(1.-2.*ga1)) &
                     +w11*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w11+w11**2)/(ga2**2+w11*(1.-2.*ga2))
                ww21m=w21*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w21+w21**2)/(ga2**2+w21*(1.-2.*ga2)) &
                     +w21*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w21+w21**2)/(ga1**2+w21*(1.-2.*ga1))
                swm=ww11m+ww21m
                w11=ww11m/swm
                w21=ww21m/swm
!
                ww12=0.5*(1.+sign(1.,v1))*(ga1/beta21) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(ga2/beta22) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
                ww12m=w12*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w12+w12**2)/(ga1**2+w12*(1.-2.*ga1)) &
                     +w12*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w12+w12**2)/(ga2**2+w12*(1.-2.*ga2))
                ww22m=w22*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w22+w22**2)/(ga2**2+w22*(1.-2.*ga2)) &
                     +w22*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w22+w22**2)/(ga1**2+w22*(1.-2.*ga1))
                swm=ww12m+ww22m
                w12=ww12m/swm
                w22=ww22m/swm
!
                ww13=0.5*(1.+sign(1.,v1))*(ga1/beta31) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(ga2/beta32) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
                ww13m=w13*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w13+w13**2)/(ga1**2+w13*(1.-2.*ga1)) &
                     +w13*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w13+w13**2)/(ga2**2+w13*(1.-2.*ga2))
                ww23m=w23*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w23+w23**2)/(ga2**2+w23*(1.-2.*ga2)) &
                     +w23*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w23+w23**2)/(ga1**2+w23*(1.-2.*ga1))
                swm=ww13m+ww23m
                w13=ww13m/swm
                w23=ww23m/swm
!
                ww14=0.5*(1.+sign(1.,v4))*(ga1/beta41) &
                     +0.5*(1.-sign(1.,v4))*(ga2/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(ga2/beta42) &
                     +0.5*(1.-sign(1.,v4))*(ga1/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
                ww14m=w14*0.5*(1.+sign(1.,v4))*(ga1+ga1**2-3.*ga1*w14+w14**2)/(ga1**2+w14*(1.-2.*ga1)) &
                     +w14*0.5*(1.-sign(1.,v4))*(ga2+ga2**2-3.*ga2*w14+w14**2)/(ga2**2+w14*(1.-2.*ga2))
                ww24m=w24*0.5*(1.+sign(1.,v4))*(ga2+ga2**2-3.*ga2*w24+w24**2)/(ga2**2+w24*(1.-2.*ga2)) &
                     +w24*0.5*(1.-sign(1.,v4))*(ga1+ga1**2-3.*ga1*w24+w24**2)/(ga1**2+w24*(1.-2.*ga1))
                swm=ww14m+ww24m
                w14=ww14m/swm
                w24=ww24m/swm
!
                ww15=0.5*(1.+sign(1.,v5))*(ga1/beta51) &
                     +0.5*(1.-sign(1.,v5))*(ga2/beta51)
                ww25=0.5*(1.+sign(1.,v5))*(ga2/beta52) &
                     +0.5*(1.-sign(1.,v5))*(ga1/beta52)
                sw=ww15+ww25
                w15=ww15/sw
                w25=ww25/sw
                ww15m=w15*0.5*(1.+sign(1.,v5))*(ga1+ga1**2-3.*ga1*w15+w15**2)/(ga1**2+w15*(1.-2.*ga1)) &
                     +w15*0.5*(1.-sign(1.,v5))*(ga2+ga2**2-3.*ga2*w15+w15**2)/(ga2**2+w15*(1.-2.*ga2))
                ww25m=w25*0.5*(1.+sign(1.,v5))*(ga2+ga2**2-3.*ga2*w25+w25**2)/(ga2**2+w25*(1.-2.*ga2)) &
                     +w25*0.5*(1.-sign(1.,v5))*(ga1+ga1**2-3.*ga1*w25+w25**2)/(ga1**2+w25*(1.-2.*ga1))
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
!        splitting du flux de l'energie
                qp5f1m=q51*0.5*(v(n-ninc,3)+abs(v(n-ninc,3))) &
                     +q52*0.5*(fxy(m-ninc)+v(n-ninc,2)*abs(v(n-ninc,3)/v(n-ninc,1))) &
                     +q53*0.5*(fyy(m-ninc)+v(n-ninc,3)*abs(v(n-ninc,3)/v(n-ninc,1))) &
                     +q54*0.5*(fyz(m-ninc)+v(n-ninc,4)*abs(v(n-ninc,3)/v(n-ninc,1) &
                     +sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1))))) &
                     +q55*0.5*(fey(m-ninc)+v(n-ninc,5)*abs(v(n-ninc,3)/v(n-ninc,1) &
                     -sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1)))))
                qm5f1m=q51*0.5*(v(n-ninc,3)-abs(v(n-ninc,3))) &
                     +q52*0.5*(fxy(m-ninc)-v(n-ninc,2)*abs(v(n-ninc,3)/v(n-ninc,1))) &
                     +q53*0.5*(fyy(m-ninc)-v(n-ninc,3)*abs(v(n-ninc,3)/v(n-ninc,1))) &
                     +q54*0.5*(fyz(m-ninc)-v(n-ninc,4)*abs(v(n-ninc,3)/v(n-ninc,1) &
                     +sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1))))) &
                     +q55*0.5*(fey(m-ninc)-v(n-ninc,5)*abs(v(n-ninc,3)/v(n-ninc,1) &
                     -sqrt(abs(gam*ps(n-ninc)/v(n-ninc,1)))))
!
                qp5f=q51*0.5*(v(n,3)+abs(v(n,3))) &
                     +q52*0.5*(fxy(m)+v(n,2)*abs(v(n,3)/v(n,1))) &
                     +q53*0.5*(fyy(m)+v(n,3)*abs(v(n,3)/v(n,1))) &
                     +q54*0.5*(fyz(m)+v(n,4)*abs(v(n,3)/v(n,1)+sqrt(abs(gam*ps(n)/v(n,1))))) &
                     +q55*0.5*(fey(m)+v(n,5)*abs(v(n,3)/v(n,1)-sqrt(abs(gam*ps(n)/v(n,1)))))
                qm5f=q51*0.5*(v(n,3)-abs(v(n,3))) &
                     +q52*0.5*(fxy(m)-v(n,2)*abs(v(n,3)/v(n,1))) &
                     +q53*0.5*(fyy(m)-v(n,3)*abs(v(n,3)/v(n,1))) &
                     +q54*0.5*(fyz(m)-v(n,4)*abs(v(n,3)/v(n,1)+sqrt(abs(gam*ps(n)/v(n,1))))) &
                     +q55*0.5*(fey(m)-v(n,5)*abs(v(n,3)/v(n,1)-sqrt(abs(gam*ps(n)/v(n,1)))))
!
                qp5f1p=q51*0.5*(v(n+ninc,3)+abs(v(n+ninc,3))) &
                     +q52*0.5*(fxy(m+ninc)+v(n+ninc,2)*abs(v(n+ninc,3)/v(n+ninc,1))) &
                     +q53*0.5*(fyy(m+ninc)+v(n+ninc,3)*abs(v(n+ninc,3)/v(n+ninc,1))) &
                     +q54*0.5*(fyz(m+ninc)+v(n+ninc,4)*abs(v(n+ninc,3)/v(n+ninc,1) &
                     +sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1))))) &
                     +q55*0.5*(fey(m+ninc)+v(n+ninc,5)*abs(v(n+ninc,3)/v(n+ninc,1) &
                     -sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1)))))
                qm5f1p=q51*0.5*(v(n+ninc,3)-abs(v(n+ninc,3))) &
                     +q52*0.5*(fxy(m+ninc)-v(n+ninc,2)*abs(v(n+ninc,3)/v(n+ninc,1))) &
                     +q53*0.5*(fyy(m+ninc)-v(n+ninc,3)*abs(v(n+ninc,3)/v(n+ninc,1))) &
                     +q54*0.5*(fyz(m+ninc)-v(n+ninc,4)*abs(v(n+ninc,3)/v(n+ninc,1) &
                     +sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1))))) &
                     +q55*0.5*(fey(m+ninc)-v(n+ninc,5)*abs(v(n+ninc,3)/v(n+ninc,1) &
                     -sqrt(abs(gam*ps(n+ninc)/v(n+ninc,1)))))
!
                qp5f2p=q51*0.5*(v(n+2*ninc,3)+abs(v(n+2*ninc,3))) &
                     +q52*0.5*(fxy(m+2*ninc)+v(n+2*ninc,2)*abs(v(n+2*ninc,3)/v(n+2*ninc,1))) &
                     +q53*0.5*(fyy(m+2*ninc)+v(n+2*ninc,3)*abs(v(n+2*ninc,3)/v(n+2*ninc,1))) &
                     +q54*0.5*(fyz(m+2*ninc)+v(n+2*ninc,4)*abs(v(n+2*ninc,3)/v(n+2*ninc,1) &
                     +sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1))))) &
                     +q55*0.5*(fey(m+2*ninc)+v(n+2*ninc,5)*abs(v(n+2*ninc,3)/v(n+2*ninc,1) &
                     -sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1)))))
                qm5f2p=q51*0.5*(v(n+2*ninc,3)-abs(v(n+2*ninc,3))) &
                     +q52*0.5*(fxy(m+2*ninc)-v(n+2*ninc,2)*abs(v(n+2*ninc,3)/v(n+2*ninc,1))) &
                     +q53*0.5*(fyy(m+2*ninc)-v(n+2*ninc,3)*abs(v(n+2*ninc,3)/v(n+2*ninc,1))) &
                     +q54*0.5*(fyz(m+2*ninc)-v(n+2*ninc,4)*abs(v(n+2*ninc,3)/v(n+2*ninc,1) &
                     +sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1))))) &
                     +q55*0.5*(fey(m+2*ninc)-v(n+2*ninc,5)*abs(v(n+2*ninc,3)/v(n+2*ninc,1) &
                     -sqrt(abs(gam*ps(n+2*ninc)/v(n+2*ninc,1)))))
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
                gp51=qp5f1m*c10 +qp5f  *c11
                gm51=qm5f  *c00 +qm5f1p*c01
                gp52=qp5f  *c00 +qp5f1p*c01
                gm52=qm5f1p*c11 +qm5f2p*c10
!        calcul des senseurs beta (au carre)
                iexp=2
!         iexp=1
                beta11=(0.5*(1.+sign(1.,v1))*(q1f-q1f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f1p-q1f)**2+eps)**iexp
                beta12=(0.5*(1.+sign(1.,v1))*(q1f1p-q1f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q1f2p-q1f1p)**2+eps)**iexp
!
                beta21=(0.5*(1.+sign(1.,v1))*(q2f-q2f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f1p-q2f)**2+eps)**iexp
                beta22=(0.5*(1.+sign(1.,v1))*(q2f1p-q2f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q2f2p-q2f1p)**2+eps)**iexp
!
                beta31=(0.5*(1.+sign(1.,v1))*(q3f-q3f1m)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f1p-q3f)**2+eps)**iexp
                beta32=(0.5*(1.+sign(1.,v1))*(q3f1p-q3f)**2 &
                     +0.5*(1.-sign(1.,v1))*(q3f2p-q3f1p)**2+eps)**iexp
!
                beta41=(0.5*(1.+sign(1.,v4))*(q4f-q4f1m)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f1p-q4f)**2+eps)**iexp
                beta42=(0.5*(1.+sign(1.,v4))*(q4f1p-q4f)**2 &
                     +0.5*(1.-sign(1.,v4))*(q4f2p-q4f1p)**2+eps)**iexp
!
                betap51=((qp5f  -qp5f1m)**2+eps)**iexp
                betam51=((qm5f1p-qm5f  )**2+eps)**iexp
                betap52=((qp5f1p-qp5f  )**2+eps)**iexp
                betam52=((qm5f2p-qm5f1p)**2+eps)**iexp
!        calculs des poids wi
                ww11=0.5*(1.+sign(1.,v1))*(ga1/beta11) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta11)
                ww21=0.5*(1.+sign(1.,v1))*(ga2/beta12) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta12)
                sw=ww11+ww21
                w11=ww11/sw
                w21=ww21/sw
                ww11m=w11*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w11+w11**2)/(ga1**2+w11*(1.-2.*ga1)) &
                     +w11*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w11+w11**2)/(ga2**2+w11*(1.-2.*ga2))
                ww21m=w21*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w21+w21**2)/(ga2**2+w21*(1.-2.*ga2)) &
                     +w21*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w21+w21**2)/(ga1**2+w21*(1.-2.*ga1))
                swm=ww11m+ww21m
                w11=ww11m/swm
                w21=ww21m/swm
!
                ww12=0.5*(1.+sign(1.,v1))*(ga1/beta21) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta21)
                ww22=0.5*(1.+sign(1.,v1))*(ga2/beta22) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta22)
                sw=ww12+ww22
                w12=ww12/sw
                w22=ww22/sw
                ww12m=w12*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w12+w12**2)/(ga1**2+w12*(1.-2.*ga1)) &
                     +w12*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w12+w12**2)/(ga2**2+w12*(1.-2.*ga2))
                ww22m=w22*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w22+w22**2)/(ga2**2+w22*(1.-2.*ga2)) &
                     +w22*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w22+w22**2)/(ga1**2+w22*(1.-2.*ga1))
                swm=ww12m+ww22m
                w12=ww12m/swm
                w22=ww22m/swm
!
                ww13=0.5*(1.+sign(1.,v1))*(ga1/beta31) &
                     +0.5*(1.-sign(1.,v1))*(ga2/beta31)
                ww23=0.5*(1.+sign(1.,v1))*(ga2/beta32) &
                     +0.5*(1.-sign(1.,v1))*(ga1/beta32)
                sw=ww13+ww23
                w13=ww13/sw
                w23=ww23/sw
                ww13m=w13*0.5*(1.+sign(1.,v1))*(ga1+ga1**2-3.*ga1*w13+w13**2)/(ga1**2+w13*(1.-2.*ga1)) &
                     +w13*0.5*(1.-sign(1.,v1))*(ga2+ga2**2-3.*ga2*w13+w13**2)/(ga2**2+w13*(1.-2.*ga2))
                ww23m=w23*0.5*(1.+sign(1.,v1))*(ga2+ga2**2-3.*ga2*w23+w23**2)/(ga2**2+w23*(1.-2.*ga2)) &
                     +w23*0.5*(1.-sign(1.,v1))*(ga1+ga1**2-3.*ga1*w23+w23**2)/(ga1**2+w23*(1.-2.*ga1))
                swm=ww13m+ww23m
                w13=ww13m/swm
                w23=ww23m/swm
!
                ww14=0.5*(1.+sign(1.,v4))*(ga1/beta41) &
                     +0.5*(1.-sign(1.,v4))*(ga2/beta41)
                ww24=0.5*(1.+sign(1.,v4))*(ga2/beta42) &
                     +0.5*(1.-sign(1.,v4))*(ga1/beta42)
                sw=ww14+ww24
                w14=ww14/sw
                w24=ww24/sw
                ww14m=w14*0.5*(1.+sign(1.,v4))*(ga1+ga1**2-3.*ga1*w14+w14**2)/(ga1**2+w14*(1.-2.*ga1)) &
                     +w14*0.5*(1.-sign(1.,v4))*(ga2+ga2**2-3.*ga2*w14+w14**2)/(ga2**2+w14*(1.-2.*ga2))
                ww24m=w24*0.5*(1.+sign(1.,v4))*(ga2+ga2**2-3.*ga2*w24+w24**2)/(ga2**2+w24*(1.-2.*ga2)) &
                     +w24*0.5*(1.-sign(1.,v4))*(ga1+ga1**2-3.*ga1*w24+w24**2)/(ga1**2+w24*(1.-2.*ga1))
                swm=ww14m+ww24m
                w14=ww14m/swm
                w24=ww24m/swm
!
                wwp15=ga1/betap51
                wwm15=ga2/betam51
                wwp25=ga2/betap52
                wwm25=ga1/betam52
                swp=wwp15+wwp25
                swm=wwm15+wwm25
                wp15=wwp15/swp
                wp25=wwp25/swp
                wm15=wwm15/swm
                wm25=wwm25/swm
                wwp15m=wp15*(ga1+ga1**2-3.*ga1*wp15+wp15**2)/(ga1**2+wp15*(1.-2.*ga1))
                wwm15m=wm15*(ga2+ga2**2-3.*ga2*wm15+wm15**2)/(ga2**2+wm15*(1.-2.*ga2))
                wwp25m=wp25*(ga2+ga2**2-3.*ga2*wp25+wp25**2)/(ga2**2+wp25*(1.-2.*ga2))
                wwm25m=wm25*(ga1+ga1**2-3.*ga1*wm25+wm25**2)/(ga1**2+wm25*(1.-2.*ga1))
                swpm=wwp15m+wwp25m
                swmm=wwm15m+wwm25m
                wp15=wwp15m/swpm
                wp25=wwp25m/swpm
                wm15=wwm15m/swmm
                wm25=wwm25m/swmm
!        calcul des flux convectifs projetes
                gc1=w11*g11+w21*g12
                gc2=w12*g21+w22*g22
                gc3=w13*g31+w23*g32
                gc4=w14*g41+w24*g42
                gpc5=wp15*gp51+wp25*gp52
                gmc5=wm15*gm51+wm25*gm52
                gc5=gpc5+gmc5
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
                gv2=(toxx(n)+toxx(n1))*sn(m1,kdir,1) &
                     +(toxy(n)+toxy(n1))*sn(m1,kdir,2)
                gv3=(toxy(n)+toxy(n1))*sn(m1,kdir,1) &
                     +(toyy(n)+toyy(n1))*sn(m1,kdir,2)
                gv5=(toxx(n )*ul+toxy(n )*vl+qcx(n ) &
                     +toxx(n1)*ur+toxy(n1)*vr+qcx(n1))*sn(m1,kdir,1) &
                     +(toxy(n )*ul+toyy(n )*vl+qcy(n ) &
                     +toxy(n1)*ur+toyy(n1)*vr+qcy(n1))*sn(m1,kdir,2)
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

!$OMP END MASTER
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
  end subroutine sch_weno3split2
end module mod_sch_weno3split2
