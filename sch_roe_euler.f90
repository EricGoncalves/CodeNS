module mod_sch_roe_euler
implicit none
contains
      subroutine sch_roe_euler( &
                 lm,ityprk, &
                 u,v,ff, &
                 equat, &
                 sn,lgsnlt, &
                 rhol,ul,vl,wl,pl,rhor,ur,vr,wr,prr, &
                 ps)
!
!***********************************************************************
!_P                          SINUMEF
!
!_DA  DATE_C : avril 2002 - Eric Goncalves / Sinumef
!
!     ACT
!_A     Schema de Roe avec extrapolation MUSCL - equations Euler
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
double precision :: am
double precision :: am2i
double precision :: ar
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
double precision :: di1
double precision :: di2
double precision :: di3
double precision :: di4
double precision :: di5
double precision :: dj1
double precision :: dj2
double precision :: dj3
double precision :: dj4
double precision :: dj5
double precision :: dk1
double precision :: dk2
double precision :: dk3
double precision :: dk4
double precision :: dk5
double precision :: dw1
double precision :: dw2
double precision :: dw3
double precision :: dw4
double precision :: dw5
double precision :: el
double precision :: er
double precision :: f1
double precision :: f2
double precision :: f3
double precision :: f4
double precision :: f5
double precision :: fex
double precision :: fey
double precision :: fez
double precision :: fxx
double precision :: fxy
double precision :: fxz
double precision :: fyy
double precision :: fyz
double precision :: fzz
double precision :: g1
double precision :: g2
double precision :: g3
double precision :: g4
double precision :: g5
double precision :: gd
double precision :: gd1
double precision :: gd2
double precision :: h1
double precision :: h2
double precision :: h3
double precision :: h4
double precision :: h5
double precision :: hl
double precision :: hm
double precision :: hr
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
double precision :: p11
double precision :: p12
double precision :: p13
double precision :: p14
double precision :: p15
double precision :: p21
double precision :: p22
double precision :: p23
double precision :: p24
double precision :: p25
double precision :: p31
double precision :: p32
double precision :: p33
double precision :: p34
double precision :: p35
double precision :: p41
double precision :: p42
double precision :: p43
double precision :: p44
double precision :: p45
double precision :: p51
double precision :: p52
double precision :: p53
double precision :: p54
double precision :: p55
double precision :: q11
double precision :: q12
double precision :: q13
double precision :: q14
double precision :: q15
double precision :: q21
double precision :: q22
double precision :: q23
double precision :: q24
double precision :: q25
double precision :: q2l
double precision :: q2r
double precision :: q31
double precision :: q32
double precision :: q33
double precision :: q34
double precision :: q35
double precision :: q41
double precision :: q42
double precision :: q43
double precision :: q44
double precision :: q45
double precision :: q51
double precision :: q52
double precision :: q53
double precision :: q54
double precision :: q55
double precision :: rhoami
double precision :: rhoiam
double precision :: rhom
double precision :: rhomi
double precision :: si0
double precision :: si1
double precision :: si2
double precision :: si3
double precision :: si4
double precision :: sj0
double precision :: sj1
double precision :: sj2
double precision :: sj3
double precision :: sj4
double precision :: sk0
double precision :: sk1
double precision :: sk2
double precision :: sk3
double precision :: sk4
double precision :: um
double precision :: v1
double precision :: v4
double precision :: v5
double precision :: vitm2
double precision :: vm
double precision :: vn
double precision :: wm
!
!-----------------------------------------------------------------------
!
      real nx,ny,nz
      character(len=7 ) :: equat
      dimension v(ip11,ip60),u(ip11,ip60),ff(ip11,ip60)
      dimension sn(lgsnlt,nind,ndir)
      dimension ps(ip11)
      dimension rhol(ip00),ul(ip00),vl(ip00),wl(ip00),pl(ip00), &
             rhor(ip00),ur(ip00),vr(ip00),wr(ip00),prr(ip00)
!
      indc(i,j,k)=n0c+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd

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
         dfxx=rhor(m)*ur(m)**2+prr(m)-pinfl+rhol(m)*ul(m)**2+pl(m)-pinfl
         dfxy=rhor(m)*ur(m)*vr(m)+rhol(m)*ul(m)*vl(m)
         dfxz=rhor(m)*ur(m)*wr(m)+rhol(m)*ul(m)*wl(m)
         dfyy=rhor(m)*vr(m)**2+prr(m)-pinfl+rhol(m)*vl(m)**2+pl(m)-pinfl
         dfyz=rhor(m)*vr(m)*wr(m)+rhol(m)*vl(m)*wl(m)
         dfzz=rhor(m)*wr(m)**2+prr(m)-pinfl+rhol(m)*wl(m)**2+pl(m)-pinfl
         dfex=(rhor(m)*er+prr(m)-pinfl)*ur(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*ul(m)
         dfey=(rhor(m)*er+prr(m)-pinfl)*vr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*vl(m)
         dfez=(rhor(m)*er+prr(m)-pinfl)*wr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*wl(m)
         f1=(rhor(m)*ur(m)+rhol(m)*ul(m))*sn(m,kdir,1) &
           +(rhor(m)*vr(m)+rhol(m)*vl(m))*sn(m,kdir,2) &
           +(rhor(m)*wr(m)+rhol(m)*wl(m))*sn(m,kdir,3)-di1
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
        fxx=v(n1,2)*(v(n1,2)/v(n1,1))+ps(n-ninc)-pinfl
        fxy=v(n1,3)*(v(n1,2)/v(n1,1))
        fxz=v(n1,4)*(v(n1,2)/v(n1,1))
        fyy=v(n1,3)*(v(n1,3)/v(n1,1))+ps(n-ninc)-pinfl
        fyz=v(n1,4)*(v(n1,3)/v(n1,1))
        fzz=v(n1,4)*(v(n1,4)/v(n1,1))+ps(n-ninc)-pinfl
        fex=(v(n1,5)+ps(n-ninc)-pinfl)*v(n1,2)/v(n1,1)
        fey=(v(n1,5)+ps(n-ninc)-pinfl)*v(n1,3)/v(n1,1)
        fez=(v(n1,5)+ps(n-ninc)-pinfl)*v(n1,4)/v(n1,1)
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
        fxx=v(n,2)*(v(n,2)/v(n,1))+ps(n)-pinfl
        fxy=v(n,3)*(v(n,2)/v(n,1))
        fxz=v(n,4)*(v(n,2)/v(n,1))
        fyy=v(n,3)*(v(n,3)/v(n,1))+ps(n)-pinfl
        fyz=v(n,4)*(v(n,3)/v(n,1))
        fzz=v(n,4)*(v(n,4)/v(n,1))+ps(n)-pinfl
        fex=(v(n,5)+ps(n)-pinfl)*v(n,2)/v(n,1)
        fey=(v(n,5)+ps(n)-pinfl)*v(n,3)/v(n,1)
        fez=(v(n,5)+ps(n)-pinfl)*v(n,4)/v(n,1)
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
         dfxx=rhor(m)*ur(m)**2+prr(m)-pinfl+rhol(m)*ul(m)**2+pl(m)-pinfl
         dfxy=rhor(m)*ur(m)*vr(m)+rhol(m)*ul(m)*vl(m)
         dfxz=rhor(m)*ur(m)*wr(m)+rhol(m)*ul(m)*wl(m)
         dfyy=rhor(m)*vr(m)**2+prr(m)-pinfl+rhol(m)*vl(m)**2+pl(m)-pinfl
         dfyz=rhor(m)*vr(m)*wr(m)+rhol(m)*vl(m)*wl(m)
         dfzz=rhor(m)*wr(m)**2+prr(m)-pinfl+rhol(m)*wl(m)**2+pl(m)-pinfl
         dfex=(rhor(m)*er+prr(m)-pinfl)*ur(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*ul(m)
         dfey=(rhor(m)*er+prr(m)-pinfl)*vr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*vl(m)
         dfez=(rhor(m)*er+prr(m)-pinfl)*wr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*wl(m)
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
        fxx=v(n1,2)*(v(n1,2)/v(n1,1))+ps(n-ninc)-pinfl
        fxy=v(n1,3)*(v(n1,2)/v(n1,1))
        fxz=v(n1,4)*(v(n1,2)/v(n1,1))
        fyy=v(n1,3)*(v(n1,3)/v(n1,1))+ps(n-ninc)-pinfl
        fyz=v(n1,4)*(v(n1,3)/v(n1,1))
        fzz=v(n1,4)*(v(n1,4)/v(n1,1))+ps(n-ninc)-pinfl
        fex=(v(n1,5)+ps(n-ninc)-pinfl)*v(n1,2)/v(n,1)
        fey=(v(n1,5)+ps(n-ninc)-pinfl)*v(n1,3)/v(n,1)
        fez=(v(n1,5)+ps(n-ninc)-pinfl)*v(n1,4)/v(n,1)
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
        fxx=v(n,2)*(v(n,2)/v(n,1))+ps(n)-pinfl
        fxy=v(n,3)*(v(n,2)/v(n,1))
        fxz=v(n,4)*(v(n,2)/v(n,1))
        fyy=v(n,3)*(v(n,3)/v(n,1))+ps(n)-pinfl
        fyz=v(n,4)*(v(n,3)/v(n,1))
        fzz=v(n,4)*(v(n,4)/v(n,1))+ps(n)-pinfl
        fex=(v(n,5)+ps(n)-pinfl)*v(n,2)/v(n,1)
        fey=(v(n,5)+ps(n)-pinfl)*v(n,3)/v(n,1)
        fez=(v(n,5)+ps(n)-pinfl)*v(n,4)/v(n,1)
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
         dfxx=rhor(m)*ur(m)**2+prr(m)-pinfl+rhol(m)*ul(m)**2+pl(m)-pinfl
         dfxy=rhor(m)*ur(m)*vr(m)+rhol(m)*ul(m)*vl(m)
         dfxz=rhor(m)*ur(m)*wr(m)+rhol(m)*ul(m)*wl(m)
         dfyy=rhor(m)*vr(m)**2+prr(m)-pinfl+rhol(m)*vl(m)**2+pl(m)-pinfl
         dfyz=rhor(m)*vr(m)*wr(m)+rhol(m)*vl(m)*wl(m)
         dfzz=rhor(m)*wr(m)**2+prr(m)-pinfl+rhol(m)*wl(m)**2+pl(m)-pinfl
         dfex=(rhor(m)*er+prr(m)-pinfl)*ur(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*ul(m)
         dfey=(rhor(m)*er+prr(m)-pinfl)*vr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*vl(m)
         dfez=(rhor(m)*er+prr(m)-pinfl)*wr(m)+ &
              (rhol(m)*el+pl(m)-pinfl)*wl(m)
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
        fxx=v(n1,2)*(v(n1,2)/v(n1,1))+ps(n-ninc)-pinfl
        fxy=v(n1,3)*(v(n1,2)/v(n1,1))
        fxz=v(n1,4)*(v(n1,2)/v(n1,1))
        fyy=v(n1,3)*(v(n1,3)/v(n1,1))+ps(n-ninc)-pinfl
        fyz=v(n1,4)*(v(n1,3)/v(n1,1))
        fzz=v(n1,4)*(v(n1,4)/v(n1,1))+ps(n-ninc)-pinfl
        fex=(v(n1,5)+ps(n-ninc)-pinfl)*v(n1,2)/v(n,1)
        fey=(v(n1,5)+ps(n-ninc)-pinfl)*v(n1,3)/v(n,1)
        fez=(v(n1,5)+ps(n-ninc)-pinfl)*v(n1,4)/v(n,1)
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
        fxx=v(n,2)*(v(n,2)/v(n,1))+ps(n)-pinfl
        fxy=v(n,3)*(v(n,2)/v(n,1))
        fxz=v(n,4)*(v(n,2)/v(n,1))
        fyy=v(n,3)*(v(n,3)/v(n,1))+ps(n)-pinfl
        fyz=v(n,4)*(v(n,3)/v(n,1))
        fzz=v(n,4)*(v(n,4)/v(n,1))+ps(n)-pinfl
        fex=(v(n,5)+ps(n)-pinfl)*v(n,2)/v(n,1)
        fey=(v(n,5)+ps(n)-pinfl)*v(n,3)/v(n,1)
        fez=(v(n,5)+ps(n)-pinfl)*v(n,4)/v(n,1)
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
end module
