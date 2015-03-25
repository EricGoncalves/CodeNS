module mod_sch_roe_prcd
  implicit none
contains
  subroutine sch_roe_prcd( &
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
!_DA  DATE_C : avril 2002 - Eric Goncalves / Sinumef
!
!     ACT
!_A     Schema de Roe avec extrapolation MUSCL.
!_A     Preconditionnement basse vitesse de Turkel (P,u,S).
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
  integer          ::       i,     i1,   i1m1,   i1p1,     i2
  integer          ::    i2m1,     id,   ind1,   ind2,isortie
  integer          ::  ityprk,      j,     j1,   j1m1,   j1p1
  integer          ::      j2,   j2m1,     jd,      k,     k1
  integer          ::    k1m1,   k1p1,     k2,   k2m1,     kd
  integer          ::    kdir, lgsnlt,     lm,      m,      n
  integer          ::     n0c,     n1,    nci,    ncj,    nck
  integer          ::     nid,   nijd,   ninc,    njd
  double precision ::                   al,                  am,                am2i,                  ar,               beta2
  double precision ::                 cnds,                dfex,                dfey,                dfez,                dfxx
  double precision ::                 dfxy,                dfxz,                dfyy,                dfyz,                dfzz
  double precision ::                  di1,                 di2,                 di3,                 di4,                 di5
  double precision ::                  dj1,                 dj2,                 dj3,                 dj4,                 dj5
  double precision ::                  dk1,                 dk2,                 dk3,                 dk4,                 dk5
  double precision ::                  dw1,                 dw2,                 dw3,                 dw4,                 dw5
  double precision ::                   el,                  er,                  f1,                  f2,                  f3
  double precision ::                   f4,                  f5,                 fex,                 fey,                 fez
  double precision ::        ff(ip11,ip60),                 fxx,                 fxy,                 fxz,                 fyy
  double precision ::                  fyz,                 fzz,                  g1,                  g2,                  g3
  double precision ::                   g4,                  g5,                  gd,                 gd1,                 gd2
  double precision ::                   h1,                  h2,                  h3,                  h4,                  h5
  double precision ::                   hl,                  hm,                  hr,                  nx,                  ny
  double precision ::                   nz,                 p11,                 p12,                 p13,                 p14
  double precision ::                  p15,                 p21,                 p22,                 p23,                 p24
  double precision ::                  p25,                 p31,                 p32,                 p33,                 p34
  double precision ::                  p35,                 p41,                 p42,                 p43,                 p44
  double precision ::                  p45,                 p51,                 p52,                 p53,                 p54
  double precision ::                  p55,            pl(ip00),           prr(ip00),            ps(ip11),                 q11
  double precision ::                  q12,                 q13,                 q14,                 q15,                 q21
  double precision ::                  q22,                 q23,                 q24,                 q25,                 q2l
  double precision ::                  q2r,                 q31,                 q32,                 q33,                 q34
  double precision ::                  q35,                 q41,                 q42,                 q43,                 q44
  double precision ::                  q45,                 q51,                 q52,                 q53,                 q54
  double precision ::                  q55,           qcx(ip12),           qcy(ip12),           qcz(ip12),                qinf
  double precision ::               rhoami,              rhoiam,          rhol(ip00),                rhom,               rhomi
  double precision ::           rhor(ip00),                  rr,                 si1,                 si2,                 si3
  double precision ::                  si4,                 si5,                 sj1,                 sj2,                 sj3
  double precision ::                  sj4,                 sj5,                 sk1,                 sk2,                 sk3
  double precision ::                  sk4,                 sk5,sn(lgsnlt,nind,ndir),                  ss,          toxx(ip12)
  double precision ::           toxy(ip12),          toxz(ip12),          toyy(ip12),          toyz(ip12),          tozz(ip12)
  double precision ::                   tt,        u(ip11,ip60),            ul(ip00),                  um,            ur(ip00)
  double precision ::         v(ip11,ip60),                  v1,                  v4,                  v5,               vitm2
  double precision ::             vl(ip00),                  vm,                  vn,            vr(ip00),            wl(ip00)
  double precision ::                   wm,            wr(ip00)
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
!


!
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
!-----calcul des densites de flux convectives +visqueuses-------------------------
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
    enddo
!
!*******************************************************************************
! calcul du flux numerique par direction suivant les etapes successives :
!    1) evaluation des variables primitives extrapolees
!    2) evaluation des matrices de passage et des valeurs
!       propres associees a la matrice jacobienne A
!       (ces quantites sont evaluees en l'etat moyen de Roe)
!    3) evaluation de la matrice |A| en l'etat moyen de Roe
!       et obtention du terme de dissipation numerique
!    4) evaluation du flux numerique
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
!        calcul des etats moyens de Roe
             gd=sqrt(rhor(m)/rhol(m))
             gd1=1./(1.+gd)
             gd2=gd*gd1
             rhom=sqrt(rhol(m)*rhor(m))
             rhomi=1./rhom
             um=gd1*ul(m)+gd2*ur(m)
             vm=gd1*vl(m)+gd2*vr(m)
             wm=gd1*wl(m)+gd2*wr(m)
             hm=gd1*hl+gd2*hr
             vitm2=0.5*(um**2+vm**2+wm**2)
             am=sqrt(gam1*(hm-vitm2))
             am2i=1./(am*am)
             vn=um*nx+vm*ny+wm*nz
             rhoiam=rhom/am
             rhoami=am2i/rhoiam
             beta2=min(max(2*vitm2*am2i,cte*qinf**2*am2i),1.)
!        valeurs propres
             v4=0.5*((1.+beta2)*vn + &
                  sqrt(((1.-beta2)*vn)**2+4.*beta2/am2i))
             v5=0.5*((1.+beta2)*vn - &
                  sqrt(((1.-beta2)*vn)**2+4.*beta2/am2i))
             rr=v4-vn*beta2
             ss=v5-vn*beta2
             tt=0.5*(v5-v4)
!        calcul des coefficients des matrices de passage a gauche Q et a droite P
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
             q41=gam1*vitm2*rhoami*ss/(beta2*tt)+vn/(tt*rhoiam)
             q42=-nx/(rhoiam*tt)-gam1*um*rhoami*ss/(beta2*tt)
             q43=-ny/(rhoiam*tt)-gam1*vm*rhoami*ss/(beta2*tt)
             q44=-nz/(rhoiam*tt)-gam1*wm*rhoami*ss/(beta2*tt)
             q45=gam1*rhoami*ss/(beta2*tt)
             q51=-gam1*vitm2*rhoami*rr/(beta2*tt)-vn/(tt*rhoiam)
             q52=nx/(rhoiam*tt)+gam1*um*rhoami*rr/(beta2*tt)
             q53=ny/(rhoiam*tt)+gam1*vm*rhoami*rr/(beta2*tt)
             q54=nz/(rhoiam*tt)+gam1*wm*rhoami*rr/(beta2*tt)
             q55=-gam1*rhoami*rr/(beta2*tt)
!
             p11=nx
             p12=ny
             p13=nz
             p14=0.5*rhoiam
             p15=0.5*rhoiam
             p21=um*nx
             p22=um*ny-rhom*nz
             p23=um*nz+rhom*ny
             p24=0.5*rhoiam*(um+nx*rr)
             p25=0.5*rhoiam*(um+nx*ss)
             p31=vm*nx+rhom*nz
             p32=vm*ny
             p33=vm*nz-rhom*nx
             p34=0.5*rhoiam*(vm+ny*rr)
             p35=0.5*rhoiam*(vm+ny*ss)
             p41=wm*nx-rhom*ny
             p42=wm*ny+rhom*nx
             p43=wm*nz
             p44=0.5*rhoiam*(wm+nz*rr)
             p45=0.5*rhoiam*(wm+nz*ss)
             p51=vitm2*nx+rhom*(vm*nz-wm*ny)
             p52=vitm2*ny+rhom*(wm*nx-um*nz)
             p53=vitm2*nz+rhom*(um*ny-vm*nx)
             p54=0.5*rhoiam*(hm+rr*vn)
             p55=0.5*rhoiam*(hm+ss*vn)
!        valeurs propres
             v1=abs(vn*cnds)
             v4=0.5*abs((1.+beta2)*vn*cnds + &
                  sqrt(((1.-beta2)*v1)**2+4.*beta2*cnds**2/am2i))
             v5=0.5*abs((1.+beta2)*vn*cnds - &
                  sqrt(((1.-beta2)*v1)**2+4.*beta2*cnds**2/am2i))
!        evaluation du terme de dissipation
             dw1=rhor(m)-rhol(m)
             dw2=rhor(m)*ur(m)-rhol(m)*ul(m)
             dw3=rhor(m)*vr(m)-rhol(m)*vl(m)
             dw4=rhor(m)*wr(m)-rhol(m)*wl(m)
             dw5=rhor(m)*er-rhol(m)*el
             di1=(p11*v1*q11+p12*v1*q21+p13*v1*q31+p14*v4*q41+p15*v5*q51)*dw1 &
                  +(p11*v1*q12+p12*v1*q22+p13*v1*q32+p14*v4*q42+p15*v5*q52)*dw2 &
                  +(p11*v1*q13+p12*v1*q23+p13*v1*q33+p14*v4*q43+p15*v5*q53)*dw3 &
                  +(p11*v1*q14+p12*v1*q24+p13*v1*q34+p14*v4*q44+p15*v5*q54)*dw4 &
                  +(p11*v1*q15+p12*v1*q25+p13*v1*q35+p14*v4*q45+p15*v5*q55)*dw5
             di2=(p21*v1*q11+p22*v1*q21+p23*v1*q31+p24*v4*q41+p25*v5*q51)*dw1 &
                  +(p21*v1*q12+p22*v1*q22+p23*v1*q32+p24*v4*q42+p25*v5*q52)*dw2 &
                  +(p21*v1*q13+p22*v1*q23+p23*v1*q33+p24*v4*q43+p25*v5*q53)*dw3 &
                  +(p21*v1*q14+p22*v1*q24+p23*v1*q34+p24*v4*q44+p25*v5*q54)*dw4 &
                  +(p21*v1*q15+p22*v1*q25+p23*v1*q35+p24*v4*q45+p25*v5*q55)*dw5
             di3=(p31*v1*q11+p32*v1*q21+p33*v1*q31+p34*v4*q41+p35*v5*q51)*dw1 &
                  +(p31*v1*q12+p32*v1*q22+p33*v1*q32+p34*v4*q42+p35*v5*q52)*dw2 &
                  +(p31*v1*q13+p32*v1*q23+p33*v1*q33+p34*v4*q43+p35*v5*q53)*dw3 &
                  +(p31*v1*q14+p32*v1*q24+p33*v1*q34+p34*v4*q44+p35*v5*q54)*dw4 &
                  +(p31*v1*q15+p32*v1*q25+p33*v1*q35+p34*v4*q45+p35*v5*q55)*dw5
             di4=(p41*v1*q11+p42*v1*q21+p43*v1*q31+p44*v4*q41+p45*v5*q51)*dw1 &
                  +(p41*v1*q12+p42*v1*q22+p43*v1*q32+p44*v4*q42+p45*v5*q52)*dw2 &
                  +(p41*v1*q13+p42*v1*q23+p43*v1*q33+p44*v4*q43+p45*v5*q53)*dw3 &
                  +(p41*v1*q14+p42*v1*q24+p43*v1*q34+p44*v4*q44+p45*v5*q54)*dw4 &
                  +(p41*v1*q15+p42*v1*q25+p43*v1*q35+p44*v4*q45+p45*v5*q55)*dw5
             di5=(p51*v1*q11+p52*v1*q21+p53*v1*q31+p54*v4*q41+p55*v5*q51)*dw1 &
                  +(p51*v1*q12+p52*v1*q22+p53*v1*q32+p54*v4*q42+p55*v5*q52)*dw2 &
                  +(p51*v1*q13+p52*v1*q23+p53*v1*q33+p54*v4*q43+p55*v5*q53)*dw3 &
                  +(p51*v1*q14+p52*v1*q24+p53*v1*q34+p54*v4*q44+p55*v5*q54)*dw4 &
                  +(p51*v1*q15+p52*v1*q25+p53*v1*q35+p54*v4*q45+p55*v5*q55)*dw5
!        calcul du flux numerique
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
!        calcul des etats moyens de Roe
             gd=sqrt(rhor(m)/rhol(m))
             gd1=1./(1.+gd)
             gd2=gd*gd1
             rhom=sqrt(rhol(m)*rhor(m))
             rhomi=1./rhom
             um=gd1*ul(m)+gd2*ur(m)
             vm=gd1*vl(m)+gd2*vr(m)
             wm=gd1*wl(m)+gd2*wr(m)
             hm=gd1*hl+gd2*hr
             vitm2=0.5*(um**2+vm**2+wm**2)
             am=sqrt(abs(gam1*(hm-vitm2)))
             am2i=1./(am*am)
             vn=um*nx+vm*ny+wm*nz
             rhoiam=rhom/am
             rhoami=am2i/rhoiam
             beta2=min(max(2*vitm2*am2i,cte*qinf**2*am2i),1.)
!        valeurs propres
             v4=0.5*((1.+beta2)*vn + &
                  sqrt(((1.-beta2)*vn)**2+4.*beta2/am2i))
             v5=0.5*((1.+beta2)*vn - &
                  sqrt(((1.-beta2)*vn)**2+4.*beta2/am2i))
             rr=v4-vn*beta2
             ss=v5-vn*beta2
             tt=0.5*(v5-v4)
!        calcul des coefficients des matrices de passage a gauche et a droite
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
             q41=gam1*vitm2*rhoami*ss/(beta2*tt)+vn/(tt*rhoiam)
             q42=-nx/(rhoiam*tt)-gam1*um*rhoami*ss/(beta2*tt)
             q43=-ny/(rhoiam*tt)-gam1*vm*rhoami*ss/(beta2*tt)
             q44=-nz/(rhoiam*tt)-gam1*wm*rhoami*ss/(beta2*tt)
             q45=gam1*rhoami*ss/(beta2*tt)
             q51=-gam1*vitm2*rhoami*rr/(beta2*tt)-vn/(tt*rhoiam)
             q52=nx/(rhoiam*tt)+gam1*um*rhoami*rr/(beta2*tt)
             q53=ny/(rhoiam*tt)+gam1*vm*rhoami*rr/(beta2*tt)
             q54=nz/(rhoiam*tt)+gam1*wm*rhoami*rr/(beta2*tt)
             q55=-gam1*rhoami*rr/(beta2*tt)
!
             p11=nx
             p12=ny
             p13=nz
             p14=0.5*rhoiam
             p15=0.5*rhoiam
             p21=um*nx
             p22=um*ny-rhom*nz
             p23=um*nz+rhom*ny
             p24=0.5*rhoiam*(um+nx*rr)
             p25=0.5*rhoiam*(um+nx*ss)
             p31=vm*nx+rhom*nz
             p32=vm*ny
             p33=vm*nz-rhom*nx
             p34=0.5*rhoiam*(vm+ny*rr)
             p35=0.5*rhoiam*(vm+ny*ss)
             p41=wm*nx-rhom*ny
             p42=wm*ny+rhom*nx
             p43=wm*nz
             p44=0.5*rhoiam*(wm+nz*rr)
             p45=0.5*rhoiam*(wm+nz*ss)
             p51=vitm2*nx+rhom*(vm*nz-wm*ny)
             p52=vitm2*ny+rhom*(wm*nx-um*nz)
             p53=vitm2*nz+rhom*(um*ny-vm*nx)
             p54=0.5*rhoiam*(hm+rr*vn)
             p55=0.5*rhoiam*(hm+ss*vn)
!        valeurs propres
             v1=abs(vn*cnds)
             v4=0.5*abs((1.+beta2)*vn*cnds + &
                  sqrt(((1.-beta2)*v1)**2+4.*beta2*cnds**2/am2i))
             v5=0.5*abs((1.+beta2)*vn*cnds - &
                  sqrt(((1.-beta2)*v1)**2+4.*beta2*cnds**2/am2i))
!        evaluation du terme de dissipation
             dw1=rhor(m)-rhol(m)
             dw2=rhor(m)*ur(m)-rhol(m)*ul(m)
             dw3=rhor(m)*vr(m)-rhol(m)*vl(m)
             dw4=rhor(m)*wr(m)-rhol(m)*wl(m)
             dw5=rhor(m)*er-rhol(m)*el
             dj1=(p11*v1*q11+p12*v1*q21+p13*v1*q31+p14*v4*q41+p15*v5*q51)*dw1 &
                  +(p11*v1*q12+p12*v1*q22+p13*v1*q32+p14*v4*q42+p15*v5*q52)*dw2 &
                  +(p11*v1*q13+p12*v1*q23+p13*v1*q33+p14*v4*q43+p15*v5*q53)*dw3 &
                  +(p11*v1*q14+p12*v1*q24+p13*v1*q34+p14*v4*q44+p15*v5*q54)*dw4 &
                  +(p11*v1*q15+p12*v1*q25+p13*v1*q35+p14*v4*q45+p15*v5*q55)*dw5
             dj2=(p21*v1*q11+p22*v1*q21+p23*v1*q31+p24*v4*q41+p25*v5*q51)*dw1 &
                  +(p21*v1*q12+p22*v1*q22+p23*v1*q32+p24*v4*q42+p25*v5*q52)*dw2 &
                  +(p21*v1*q13+p22*v1*q23+p23*v1*q33+p24*v4*q43+p25*v5*q53)*dw3 &
                  +(p21*v1*q14+p22*v1*q24+p23*v1*q34+p24*v4*q44+p25*v5*q54)*dw4 &
                  +(p21*v1*q15+p22*v1*q25+p23*v1*q35+p24*v4*q45+p25*v5*q55)*dw5
             dj3=(p31*v1*q11+p32*v1*q21+p33*v1*q31+p34*v4*q41+p35*v5*q51)*dw1 &
                  +(p31*v1*q12+p32*v1*q22+p33*v1*q32+p34*v4*q42+p35*v5*q52)*dw2 &
                  +(p31*v1*q13+p32*v1*q23+p33*v1*q33+p34*v4*q43+p35*v5*q53)*dw3 &
                  +(p31*v1*q14+p32*v1*q24+p33*v1*q34+p34*v4*q44+p35*v5*q54)*dw4 &
                  +(p31*v1*q15+p32*v1*q25+p33*v1*q35+p34*v4*q45+p35*v5*q55)*dw5
             dj4=(p41*v1*q11+p42*v1*q21+p43*v1*q31+p44*v4*q41+p45*v5*q51)*dw1 &
                  +(p41*v1*q12+p42*v1*q22+p43*v1*q32+p44*v4*q42+p45*v5*q52)*dw2 &
                  +(p41*v1*q13+p42*v1*q23+p43*v1*q33+p44*v4*q43+p45*v5*q53)*dw3 &
                  +(p41*v1*q14+p42*v1*q24+p43*v1*q34+p44*v4*q44+p45*v5*q54)*dw4 &
                  +(p41*v1*q15+p42*v1*q25+p43*v1*q35+p44*v4*q45+p45*v5*q55)*dw5
             dj5=(p51*v1*q11+p52*v1*q21+p53*v1*q31+p54*v4*q41+p55*v5*q51)*dw1 &
                  +(p51*v1*q12+p52*v1*q22+p53*v1*q32+p54*v4*q42+p55*v5*q52)*dw2 &
                  +(p51*v1*q13+p52*v1*q23+p53*v1*q33+p54*v4*q43+p55*v5*q53)*dw3 &
                  +(p51*v1*q14+p52*v1*q24+p53*v1*q34+p54*v4*q44+p55*v5*q54)*dw4 &
                  +(p51*v1*q15+p52*v1*q25+p53*v1*q35+p54*v4*q45+p55*v5*q55)*dw5
!        calcul du flux numerique
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
!c------direction k-------------------------------------------------------
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
!        calcul des etats moyens de Roe
                gd=sqrt(rhor(m)/rhol(m))
                gd1=1./(1.+gd)
                gd2=gd*gd1
                rhom=sqrt(rhol(m)*rhor(m))
                rhomi=1./rhom
                um=gd1*ul(m)+gd2*ur(m)
                vm=gd1*vl(m)+gd2*vr(m)
                wm=gd1*wl(m)+gd2*wr(m)
                hm=gd1*hl+gd2*hr
                vitm2=0.5*(um**2+vm**2+wm**2)
                am=sqrt(abs(gam1*(hm-vitm2)))
                am2i=1./(am*am)
                vn=um*nx+vm*ny+wm*nz
                rhoiam=rhom/am
                rhoami=am2i/rhoiam
                beta2=min(max(2*vitm2*am2i,cte*qinf**2*am2i),1.)
!        valeurs propres
                v4=0.5*((1.+beta2)*vn + &
                     sqrt(((1.-beta2)*vn)**2+4.*beta2/am2i))
                v5=0.5*((1.+beta2)*vn - &
                     sqrt(((1.-beta2)*vn)**2+4.*beta2/am2i))
                rr=v4-vn*beta2
                ss=v5-vn*beta2
                tt=0.5*(v5-v4)
!        calcul des coefficients des matrices de passage a gauche et a droite
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
                q41=gam1*vitm2*rhoami*ss/(beta2*tt)+vn/(tt*rhoiam)
                q42=-nx/(rhoiam*tt)-gam1*um*rhoami*ss/(beta2*tt)
                q43=-ny/(rhoiam*tt)-gam1*vm*rhoami*ss/(beta2*tt)
                q44=-nz/(rhoiam*tt)-gam1*wm*rhoami*ss/(beta2*tt)
                q45=gam1*rhoami*ss/(beta2*tt)
                q51=-gam1*vitm2*rhoami*rr/(beta2*tt)-vn/(tt*rhoiam)
                q52=nx/(rhoiam*tt)+gam1*um*rhoami*rr/(beta2*tt)
                q53=ny/(rhoiam*tt)+gam1*vm*rhoami*rr/(beta2*tt)
                q54=nz/(rhoiam*tt)+gam1*wm*rhoami*rr/(beta2*tt)
                q55=-gam1*rhoami*rr/(beta2*tt)
!
                p11=nx
                p12=ny
                p13=nz
                p14=0.5*rhoiam
                p15=0.5*rhoiam
                p21=um*nx
                p22=um*ny-rhom*nz
                p23=um*nz+rhom*ny
                p24=0.5*rhoiam*(um+nx*rr)
                p25=0.5*rhoiam*(um+nx*ss)
                p31=vm*nx+rhom*nz
                p32=vm*ny
                p33=vm*nz-rhom*nx
                p34=0.5*rhoiam*(vm+ny*rr)
                p35=0.5*rhoiam*(vm+ny*ss)
                p41=wm*nx-rhom*ny
                p42=wm*ny+rhom*nx
                p43=wm*nz
                p44=0.5*rhoiam*(wm+nz*rr)
                p45=0.5*rhoiam*(wm+nz*ss)
                p51=vitm2*nx+rhom*(vm*nz-wm*ny)
                p52=vitm2*ny+rhom*(wm*nx-um*nz)
                p53=vitm2*nz+rhom*(um*ny-vm*nx)
                p54=0.5*rhoiam*(hm+rr*vn)
                p55=0.5*rhoiam*(hm+ss*vn)
!        valeurs propres
                v1=abs(vn*cnds)
                v4=0.5*abs((1.+beta2)*vn*cnds + &
                     sqrt(((1.-beta2)*v1)**2+4.*beta2*cnds**2/am2i))
                v5=0.5*abs((1.+beta2)*vn*cnds - &
                     sqrt(((1.-beta2)*v1)**2+4.*beta2*cnds**2/am2i))
!        evaluation du terme de dissipation
                dw1=rhor(m)-rhol(m)
                dw2=rhor(m)*ur(m)-rhol(m)*ul(m)
                dw3=rhor(m)*vr(m)-rhol(m)*vl(m)
                dw4=rhor(m)*wr(m)-rhol(m)*wl(m)
                dw5=rhor(m)*er-rhol(m)*el
                dk1=(p11*v1*q11+p12*v1*q21+p13*v1*q31+p14*v4*q41+p15*v5*q51)*dw1 &
                     +(p11*v1*q12+p12*v1*q22+p13*v1*q32+p14*v4*q42+p15*v5*q52)*dw2 &
                     +(p11*v1*q13+p12*v1*q23+p13*v1*q33+p14*v4*q43+p15*v5*q53)*dw3 &
                     +(p11*v1*q14+p12*v1*q24+p13*v1*q34+p14*v4*q44+p15*v5*q54)*dw4 &
                     +(p11*v1*q15+p12*v1*q25+p13*v1*q35+p14*v4*q45+p15*v5*q55)*dw5
                dk2=(p21*v1*q11+p22*v1*q21+p23*v1*q31+p24*v4*q41+p25*v5*q51)*dw1 &
                     +(p21*v1*q12+p22*v1*q22+p23*v1*q32+p24*v4*q42+p25*v5*q52)*dw2 &
                     +(p21*v1*q13+p22*v1*q23+p23*v1*q33+p24*v4*q43+p25*v5*q53)*dw3 &
                     +(p21*v1*q14+p22*v1*q24+p23*v1*q34+p24*v4*q44+p25*v5*q54)*dw4 &
                     +(p21*v1*q15+p22*v1*q25+p23*v1*q35+p24*v4*q45+p25*v5*q55)*dw5
                dk3=(p31*v1*q11+p32*v1*q21+p33*v1*q31+p34*v4*q41+p35*v5*q51)*dw1 &
                     +(p31*v1*q12+p32*v1*q22+p33*v1*q32+p34*v4*q42+p35*v5*q52)*dw2 &
                     +(p31*v1*q13+p32*v1*q23+p33*v1*q33+p34*v4*q43+p35*v5*q53)*dw3 &
                     +(p31*v1*q14+p32*v1*q24+p33*v1*q34+p34*v4*q44+p35*v5*q54)*dw4 &
                     +(p31*v1*q15+p32*v1*q25+p33*v1*q35+p34*v4*q45+p35*v5*q55)*dw5
                dk4=(p41*v1*q11+p42*v1*q21+p43*v1*q31+p44*v4*q41+p45*v5*q51)*dw1 &
                     +(p41*v1*q12+p42*v1*q22+p43*v1*q32+p44*v4*q42+p45*v5*q52)*dw2 &
                     +(p41*v1*q13+p42*v1*q23+p43*v1*q33+p44*v4*q43+p45*v5*q53)*dw3 &
                     +(p41*v1*q14+p42*v1*q24+p43*v1*q34+p44*v4*q44+p45*v5*q54)*dw4 &
                     +(p41*v1*q15+p42*v1*q25+p43*v1*q35+p44*v4*q45+p45*v5*q55)*dw5
                dk5=(p51*v1*q11+p52*v1*q21+p53*v1*q31+p54*v4*q41+p55*v5*q51)*dw1 &
                     +(p51*v1*q12+p52*v1*q22+p53*v1*q32+p54*v4*q42+p55*v5*q52)*dw2 &
                     +(p51*v1*q13+p52*v1*q23+p53*v1*q33+p54*v4*q43+p55*v5*q53)*dw3 &
                     +(p51*v1*q14+p52*v1*q24+p53*v1*q34+p54*v4*q44+p55*v5*q54)*dw4 &
                     +(p51*v1*q15+p52*v1*q25+p53*v1*q35+p54*v4*q45+p55*v5*q55)*dw5
!        calcul du flux numerique
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
       write(6,'("===>sch_roe_prcd: ecriture increment expli")')
       k=1
       i=13
       do j=j1,j2m1
          n=indc(i,j,k)
          m=n-n0c
          write(6,'(i4,i6,4(1pe12.4))') &
               j,n,u(n,1),u(n,2),u(n,3),u(n,5)
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
  end subroutine sch_roe_prcd
end module mod_sch_roe_prcd
