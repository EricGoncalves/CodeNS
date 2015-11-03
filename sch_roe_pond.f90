module mod_sch_roe_pond
  implicit none
contains
  subroutine sch_roe_pond( &
       lm,ityprk, &
       u,v,ff, &
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
!_DA  DATE_C : mars 2002 - Eric Goncalves / Sinumef
!
!     ACT
!_A    Calcul des bilans de flux physiques pour chaque cellule
!_A    Schema de Roe - ordre 1 et superieur avec extrapolation MUSCL.
!_A    Ponderation pour maillage irregulier.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    use schemanum
    implicit none
    integer          ::       i,     i1,   i1m1,   i1p1,     i2
    integer          ::    i2m1,     id,   ind1,   ind2,isortie
    integer          ::  ityprk,      j,     j1,   j1m1,   j1p1
    integer          ::      j2,   j2m1,     jd,      k,     k1
    integer          ::    k1m1,   k1p1,     k2,   k2m1,     kd
    integer          ::    kdir, lgsnlt,     lm,      m,      n
    integer          ::     n0c,     n1,    nci,    ncj,    nck
    integer          ::     nid,   nijd,   ninc,    njd
    double precision ::     al,    am,  am2i,    ar,   cal
    double precision ::    car,   cbl,   cbr, cmui1, cmui2
    double precision ::  cmuj1, cmuj2, cmuk1, cmuk2,  cnds
    double precision ::    cvi,   cvj,   cvk,  dfex,  dfey
    double precision ::   dfez,  dfxx,  dfxy,  dfxz,  dfyy
    double precision ::   dfyz,  dfzz,   di1,   di2,   di3
    double precision ::    di4,   di5,   dj1,   dj2,   dj3
    double precision ::    dj4,   dj5,   dk1,   dk2,   dk3
    double precision ::    dk4,   dk5,   dw1,   dw2,   dw3
    double precision ::    dw4,   dw5,    el,    er,    f1
    double precision ::     f2,    f3,    f4,    f5,   fex
    double precision ::    fey,   fez,    ff,   fxx,   fxy
    double precision ::    fxz,   fyy,   fyz,   fzz,    g1
    double precision ::     g2,    g3,    g4,    g5,    gd
    double precision ::    gd1,   gd2,    h1,    h2,    h3
    double precision ::     h4,    h5,    hl,    hm,    hr
    double precision ::     nx,    ny,    nz,   p11,   p12
    double precision ::    p13,   p14,   p15,   p21,   p22
    double precision ::    p23,   p24,   p25,   p31,   p32
    double precision ::    p33,   p34,   p35,   p41,   p42
    double precision ::    p43,   p44,   p45,   p51,   p52
    double precision ::    p53,   p54,   p55,    pl,   prr
    double precision ::     ps,   q11,   q12,   q13,   q14
    double precision ::    q15,   q21,   q22,   q23,   q24
    double precision ::    q25,   q2l,   q2r,   q31,   q32
    double precision ::    q33,   q34,   q35,   q41,   q42
    double precision ::    q43,   q44,   q45,   q51,   q52
    double precision ::    q53,   q54,   q55,   qcx,   qcy
    double precision ::    qcz,rhoami,rhoiam,  rhol,  rhom
    double precision ::  rhomi,  rhor,   si0,   si1,   si2
    double precision ::    si3,   si4,   sj0,   sj1,   sj2
    double precision ::    sj3,   sj4,   sk0,   sk1,   sk2
    double precision ::    sk3,   sk4,    sn,  toxx,  toxy
    double precision ::   toxz,  toyy,  toyz,  tozz,     u
    double precision ::     ul,    um,    ur,     v,    v1
    double precision ::     v4,    v5, vitm2,    vl,    vm
    double precision ::     vn,    vr,    wl,    wm,    wr
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
    dimension v(ip11,ip60),u(ip11,ip60),ff(ip11,ip60)
    dimension toxx(ip12),toxy(ip12),toxz(ip12), &
         toyy(ip12),toyz(ip12),tozz(ip12), &
         qcx (ip12),qcy (ip12),qcz (ip12)
    dimension sn(lgsnlt,nind,ndir)
    dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
              cmuk1(ip21),cmuk2(ip21),cvi(ip21),cvj(ip21),cvk(ip21)
    dimension ps(ip11)
    dimension rhol(ip00),ul(ip00),vl(ip00),wl(ip00),pl(ip00), &
              rhor(ip00),ur(ip00),vr(ip00),wr(ip00),prr(ip00)
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
    isortie=0
!
!-----calcul des densites de flux convectives +visqueuses--------------------
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
!------direction i----------------------------------------------
!
    kdir=1
    ninc=nci
!
!-----definition des variables extrapolees------------------------
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1p1,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
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
             vitm2=0.5*(um*um+vm*vm+wm*wm)
             am=sqrt(gam1*(hm-vitm2))
             am2i=1./(am*am)
             vn=um*nx+vm*ny+wm*nz
             rhoiam=rhom/am
             rhoami=am2i/rhoiam
!        valeurs propres
             v1=abs(vn*cnds)
             v4=abs((vn+am)*cnds)
             v5=abs((vn-am)*cnds)
!        calcul des coefficients des matrices de passage à gauche et a droite
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
             dfxx=cmui1(n)*(rhor(m)*ur(m)**2 + prr(m)-pinfl - toxx(n)) &
                 +cmui2(n)*(rhol(m)*ul(m)**2 + pl(m)-pinfl - toxx(n-ninc))
             dfxy=cmui1(n)*(rhor(m)*ur(m)*vr(m) - toxy(n)) &
                 +cmui2(n)*(rhol(m)*ul(m)*vl(m) - toxy(n-ninc))
             dfxz=cmui1(n)*(rhor(m)*ur(m)*wr(m) - toxz(n)) &
                 +cmui2(n)*(rhol(m)*ul(m)*wl(m) - toxz(n-ninc))
             dfyy=cmui1(n)*(rhor(m)*vr(m)**2 + prr(m)-pinfl - toyy(n)) &
                 +cmui2(n)*(rhol(m)*vl(m)**2 + pl(m)-pinfl - toyy(n-ninc))
             dfyz=cmui1(n)*(rhor(m)*vr(m)*wr(m) - toyz(n)) &
                 +cmui2(n)*(rhol(m)*vl(m)*wl(m) - toyz(n-ninc))
             dfzz=cmui1(n)*(rhor(m)*wr(m)**2 + prr(m)-pinfl - tozz(n)) &
                 +cmui2(n)*(rhol(m)*wl(m)**2 + pl(m)-pinfl - tozz(n-ninc))
             dfex=cmui1(n)*((rhor(m)*er+prr(m)-pinfl)*ur(m) - qcx(n) &
                  -(toxx(n     )*ur(m)+toxy(n     )*vr(m)+toxz(n     )*wr(m))) &
                 +cmui2(n)*((rhol(m)*el+pl(m)-pinfl)*ul(m) - qcx(n-ninc) &
                  -(toxx(n-ninc)*ul(m)+toxy(n-ninc)*vl(m)+toxz(n-ninc)*wl(m)))
             dfey=cmui1(n)*((rhor(m)*er+prr(m)-pinfl)*vr(m) - qcy(n) &
                  -(toxy(n     )*ur(m)+toyy(n     )*vr(m)+toyz(n     )*wr(m))) &
                 +cmui2(n)*((rhol(m)*el+pl(m)-pinfl)*vl(m) - qcy(n-ninc) &
                  -(toxy(n-ninc)*ul(m)+toyy(n-ninc)*vl(m)+toyz(n-ninc)*wl(m)))
             dfez=cmui1(n)*((rhor(m)*er+prr(m)-pinfl)*wr(m) - qcz(n) &
                  -(toxz(n     )*ur(m)+toyz(n     )*vr(m)+tozz(n     )*wr(m))) &
                 +cmui2(n)*((rhol(m)*el+pl(m)-pinfl)*wl(m) - qcz(n-ninc) &
                  -(toxz(n-ninc)*ul(m)+toyz(n-ninc)*vl(m)+tozz(n-ninc)*wl(m)))
             f1=(cmui1(n)*rhor(m)*ur(m)+cmui2(n)*rhol(m)*ul(m))*sn(m,kdir,1) &
               +(cmui1(n)*rhor(m)*vr(m)+cmui2(n)*rhol(m)*vl(m))*sn(m,kdir,2) &
               +(cmui1(n)*rhor(m)*wr(m)+cmui2(n)*rhol(m)*wl(m))*sn(m,kdir,3) &
                  -di1
             f2=dfxx*sn(m,kdir,1)+dfxy*sn(m,kdir,2)+dfxz*sn(m,kdir,3)-di2
             f3=dfxy*sn(m,kdir,1)+dfyy*sn(m,kdir,2)+dfyz*sn(m,kdir,3)-di3
             f4=dfxz*sn(m,kdir,1)+dfyz*sn(m,kdir,2)+dfzz*sn(m,kdir,3)-di4
             f5=dfex*sn(m,kdir,1)+dfey*sn(m,kdir,2)+dfez*sn(m,kdir,3)-di5
!
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
             vitm2=0.5*(um*um+vm*vm+wm*wm)
             am=sqrt(gam1*(hm-vitm2))
             am2i=1./(am*am)
             vn=um*nx+vm*ny+wm*nz
             rhoiam=rhom/am
             rhoami=am2i/rhoiam
!        valeurs propres
             v1=abs(vn*cnds)
             v4=abs((vn+am)*cnds)
             v5=abs((vn-am)*cnds)
!        calcul des coefficients des matrices de passage à gauche et a droite
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
             dfxx=cmuj1(n)*(rhor(m)*ur(m)**2 + prr(m)-pinfl - toxx(n)) &
                 +cmuj2(n)*(rhol(m)*ul(m)**2 + pl(m)-pinfl - toxx(n-ninc))
             dfxy=cmuj1(n)*(rhor(m)*ur(m)*vr(m) - toxy(n)) &
                 +cmuj2(n)*(rhol(m)*ul(m)*vl(m) - toxy(n-ninc))
             dfxz=cmuj1(n)*(rhor(m)*ur(m)*wr(m) - toxz(n)) &
                 +cmuj2(n)*(rhol(m)*ul(m)*wl(m) - toxz(n-ninc))
             dfyy=cmuj1(n)*(rhor(m)*vr(m)**2 + prr(m)-pinfl - toyy(n)) &
                 +cmuj2(n)*(rhol(m)*vl(m)**2 + pl(m)-pinfl - toyy(n-ninc))
             dfyz=cmuj1(n)*(rhor(m)*vr(m)*wr(m) - toyz(n)) &
                 +cmuj2(n)*(rhol(m)*vl(m)*wl(m) - toyz(n-ninc))
             dfzz=cmuj1(n)*(rhor(m)*wr(m)**2 + prr(m)-pinfl - tozz(n)) &
                 +cmuj2(n)*(rhol(m)*wl(m)**2 + pl(m)-pinfl - tozz(n-ninc))
             dfex=cmuj1(n)*((rhor(m)*er+prr(m)-pinfl)*ur(m) - qcx(n) &
                  -(toxx(n     )*ur(m)+toxy(n     )*vr(m)+toxz(n     )*wr(m))) &
                 +cmuj2(n)*((rhol(m)*el+pl(m)-pinfl)*ul(m) - qcx(n-ninc) &
                  -(toxx(n-ninc)*ul(m)+toxy(n-ninc)*vl(m)+toxz(n-ninc)*wl(m)))
             dfey=cmuj1(n)*((rhor(m)*er+prr(m)-pinfl)*vr(m) - qcy(n) &
                  -(toxy(n     )*ur(m)+toyy(n     )*vr(m)+toyz(n     )*wr(m))) &
                 +cmuj2(n)*((rhol(m)*el+pl(m)-pinfl)*vl(m) - qcy(n-ninc) &
                  -(toxy(n-ninc)*ul(m)+toyy(n-ninc)*vl(m)+toyz(n-ninc)*wl(m)))
             dfez=cmuj1(n)*((rhor(m)*er+prr(m)-pinfl)*wr(m) - qcz(n) &
                  -(toxz(n     )*ur(m)+toyz(n     )*vr(m)+tozz(n     )*wr(m))) &
                 +cmuj2(n)*((rhol(m)*el+pl(m)-pinfl)*wl(m) - qcz(n-ninc) &
                  -(toxz(n-ninc)*ul(m)+toyz(n-ninc)*vl(m)+tozz(n-ninc)*wl(m)))
             g1=(cmuj1(n)*rhor(m)*ur(m)+cmuj2(n)*rhol(m)*ul(m))*sn(m,kdir,1) &
               +(cmuj1(n)*rhor(m)*vr(m)+cmuj2(n)*rhol(m)*vl(m))*sn(m,kdir,2) &
               +(cmuj1(n)*rhor(m)*wr(m)+cmuj2(n)*rhol(m)*wl(m))*sn(m,kdir,3) &
                  -dj1
             g2=dfxx*sn(m,kdir,1)+dfxy*sn(m,kdir,2)+dfxz*sn(m,kdir,3)-dj2
             g3=dfxy*sn(m,kdir,1)+dfyy*sn(m,kdir,2)+dfyz*sn(m,kdir,3)-dj3
             g4=dfxz*sn(m,kdir,1)+dfyz*sn(m,kdir,2)+dfzz*sn(m,kdir,3)-dj4
             g5=dfex*sn(m,kdir,1)+dfey*sn(m,kdir,2)+dfez*sn(m,kdir,3)-dj5
!
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
!
       do k=k1p1,k2m1
          do j=j1,j2m1
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
                am=sqrt(gam1*(hm-vitm2))
                am2i=1./(am*am)
                vn=um*nx+vm*ny+wm*nz
                rhoiam=rhom/am
                rhoami=am2i/rhoiam
!        valeurs propres
                v1=abs(vn*cnds)
                v4=abs(vn*cnds+am*cnds)
                v5=abs(vn*cnds-am*cnds)
!        calcul des coefficients des matrices de passage à gauche et a droite
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
                dfxx=cmuk1(n)*(rhor(m)*ur(m)**2 + prr(m)-pinfl - toxx(n)) &
                    +cmuk2(n)*(rhol(m)*ul(m)**2 + pl(m)-pinfl - toxx(n-ninc))
                dfxy=cmuk1(n)*(rhor(m)*ur(m)*vr(m) - toxy(n)) &
                    +cmuk2(n)*(rhol(m)*ul(m)*vl(m) - toxy(n-ninc))
                dfxz=cmuk1(n)*(rhor(m)*ur(m)*wr(m) - toxz(n)) &
                    +cmuk2(n)*(rhol(m)*ul(m)*wl(m) - toxz(n-ninc))
                dfyy=cmuk1(n)*(rhor(m)*vr(m)**2 + prr(m)-pinfl - toyy(n)) &
                    +cmuk2(n)*(rhol(m)*vl(m)**2 + pl(m)-pinfl - toyy(n-ninc))
                dfyz=cmuk1(n)*(rhor(m)*vr(m)*wr(m) - toyz(n)) &
                    +cmuk2(n)*(rhol(m)*vl(m)*wl(m) - toyz(n-ninc))
                dfzz=cmuk1(n)*(rhor(m)*wr(m)**2 + prr(m)-pinfl - tozz(n)) &
                    +cmuk2(n)*(rhol(m)*wl(m)**2 + pl(m)-pinfl - tozz(n-ninc))
                dfex=cmuk1(n)*((rhor(m)*er+prr(m)-pinfl)*ur(m) - qcx(n) &
                     -(toxx(n     )*ur(m)+toxy(n     )*vr(m)+toxz(n     )*wr(m))) &
                    +cmuk2(n)*((rhol(m)*el+pl(m)-pinfl)*ul(m) - qcx(n-ninc) &
                     -(toxx(n-ninc)*ul(m)+toxy(n-ninc)*vl(m)+toxz(n-ninc)*wl(m)))
                dfey=cmuk1(n)*((rhor(m)*er+prr(m)-pinfl)*vr(m) - qcy(n) &
                     -(toxy(n     )*ur(m)+toyy(n     )*vr(m)+toyz(n     )*wr(m))) &
                    +cmuk2(n)*((rhol(m)*el+pl(m)-pinfl)*vl(m) - qcy(n-ninc) &
                     -(toxy(n-ninc)*ul(m)+toyy(n-ninc)*vl(m)+toyz(n-ninc)*wl(m)))
                dfez=cmuk1(n)*((rhor(m)*er+prr(m)-pinfl)*wr(m) - qcz(n) &
                     -(toxz(n     )*ur(m)+toyz(n     )*vr(m)+tozz(n     )*wr(m))) &
                    +cmuk2(n)*((rhol(m)*el+pl(m)-pinfl)*wl(m) - qcz(n-ninc) &
                     -(toxz(n-ninc)*ul(m)+toyz(n-ninc)*vl(m)+tozz(n-ninc)*wl(m)))
                h1=(cmuk1(n)*rhor(m)*ur(m)+cmuk2(n)*rhol(m)*ul(m))*sn(m,kdir,1) &
                  +(cmuk1(n)*rhor(m)*vr(m)+cmuk2(n)*rhol(m)*vl(m))*sn(m,kdir,2) &
                  +(cmuk1(n)*rhor(m)*wr(m)+cmuk2(n)*rhol(m)*wl(m))*sn(m,kdir,3) &
                     -dk1
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
    if(isortie.eq.1) then
       write(6,'("===>sch_roe_pond: ecriture increment expli")')
       k=1
!      i=13
       i=160
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
  end subroutine sch_roe_pond
end module mod_sch_roe_pond
