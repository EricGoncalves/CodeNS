module mod_dissip_jameson_prcd2
  implicit none
contains
  subroutine dissip_jameson_prcd2( &
       lm,v,d, &
       equat, &
       sn,lgsnlt, &
       snup,ps,cson)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 2006 - Eric Goncalves / LEGI
!
!     ACT
!       Calcul de la dissipation scalaire du schema de Jameson.
!       Preconditionnement basse vitesse (variables P,u,e).
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
    integer          ::      i,    i1,  i1m1,  i1p1,    i2
    integer          ::   i2m1,    id,  ind1,  ind2,     j
    integer          ::     j1,  j1m1,  j1p1,    j2,  j2m1
    integer          ::     jd,     k,    k1,  k1m1,  k1p1
    integer          ::     k2,  k2m1,    kd,  kdir,lgsnlt
    integer          ::     lm,     m,     n,   n0c,   nci
    integer          ::    ncj,   nck,   nid,  nijd,  ninc
    integer          ::    njd
    double precision ::                   a2,               beta2,                 ck2,                 ck4,                cnds
    double precision ::           cson(ip11),        d(ip11,ip60),                  d1,                  d2,                  d3
    double precision ::                 d3w1,                d3w2,                d3w3,                d3w4,                d3w5
    double precision ::                   d4,                  d5,                 dd1,                 dd2,                 dd3
    double precision ::                  dd4,                 dd5,                 ds1,                 ds2,                 ds3
    double precision ::                  ds4,                 ds5,                 dw1,                 dw2,                 dw3
    double precision ::                  dw4,                 dw5,                eps2,                eps4,                  gd
    double precision ::                   ge,                 get,               prec1,               prec3,                pres
    double precision ::             ps(ip11),                  q2,                qinf,                  qq,                 rho
    double precision ::                   rl,sn(lgsnlt,nind,ndir),          snup(ip00),                  uu,        v(ip11,ip60)
    double precision ::                   vn,                  vv,                  ww
    double precision,allocatable :: temp(:)
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
!$OMP MASTER
!
    allocate(temp(ip11))


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
!c    coefficient de dissipation artificielle champ moyen
    ck2=ki2(lm)
    ck4=ki4(lm)
!
    qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
!
!-----calcul des densites de flux visqueuses--------------------------
!
    if(equat(3:5).eq.'2dk') then
       ind1 = indc(i1m1,j1m1,k1  )
       ind2 = indc(i2  ,j2  ,k2m1)
    elseif(equat(3:4).eq.'3d') then
       ind1 = indc(i1m1,j1m1,k1m1)
       ind2 = indc(i2  ,j2  ,k2  )
    endif
!$OMP SIMD
    do n=ind1,ind2
       d(n,1)=0.
       d(n,2)=0.
       d(n,3)=0.
       d(n,4)=0.
       d(n,5)=0.
       ps(n)=gam1*(v(n,5)-pinfl- &
            0.5*(v(n,2)**2+v(n,3)**2+v(n,4)**2)/v(n,1))
    enddo
!
!*********************************************************************
!      calcul du flux de dissipation artificielle par direction
!*********************************************************************
!
!------direction i----------------------------------------------
!
    kdir=1
    ninc=nci
!
!c---senseur dissipation artificielle
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             snup(m)=abs(ps(n+ninc)-2.*ps(n)+ps(n-ninc)) &
                  /(ps(n+ninc)+2.*ps(n)+ps(n-ninc))
          enddo
       enddo
    enddo
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1p1,j,k)
          ind2 = indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             eps2=ck2*max(snup(m-ninc),snup(m))
             eps4=max(0.,ck4-eps2)
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
!         a2=0.25*(cson(n)+cson(n-ninc))**2
             a2=cson(n)*cson(n-ninc)
             q2=uu**2+vv**2+ww**2
             beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
             rl=0.5*((1.+beta2)*abs(vn) + &
                  sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
             get=0.5*(v(n,5)/v(n,1)+v(n-ninc,5)/v(n-ninc,1)) !energietotale
             ge=get-0.5*q2
             gd=(1./beta2-1.)/ge  !matrice de preconditionnement Pc^-1
!
             dw1=v(n,1)-v(n-ninc,1)
             dw2=v(n,2)-v(n-ninc,2)
             dw3=v(n,3)-v(n-ninc,3)
             dw4=v(n,4)-v(n-ninc,4)
             dw5=v(n,5)-v(n-ninc,5)
             d3w1=v(n+ninc,1)-3.*v(n,1)+3.*v(n-ninc,1)-v(n-2*ninc,1)
             d3w2=v(n+ninc,2)-3.*v(n,2)+3.*v(n-ninc,2)-v(n-2*ninc,2)
             d3w3=v(n+ninc,3)-3.*v(n,3)+3.*v(n-ninc,3)-v(n-2*ninc,3)
             d3w4=v(n+ninc,4)-3.*v(n,4)+3.*v(n-ninc,4)-v(n-2*ninc,4)
             d3w5=v(n+ninc,5)-3.*v(n,5)+3.*v(n-ninc,5)-v(n-2*ninc,5)
!
             prec1=gd*(0.5*q2*dw1-uu*dw2-vv*dw3-ww*dw4+dw5)
             d1= dw1 +    prec1
             d2= dw2 + uu*prec1
             d3= dw3 + vv*prec1
             d4= dw4 + ww*prec1
             d5= dw5 +get*prec1
             prec3=gd*(0.5*q2*d3w1-uu*d3w2-vv*d3w3-ww*d3w4+d3w5)
             dd1= d3w1 +    prec3
             dd2= d3w2 + uu*prec3
             dd3= d3w3 + vv*prec3
             dd4= d3w4 + ww*prec3
             dd5= d3w5 +get*prec3
!
             ds1=eps2*rl*d1 - eps4*rl*dd1
             ds2=eps2*rl*d2 - eps4*rl*dd2
             ds3=eps2*rl*d3 - eps4*rl*dd3
             ds4=eps2*rl*d4 - eps4*rl*dd4
             ds5=eps2*rl*d5 - eps4*rl*dd5
!
             d(n,1)=d(n,1)-ds1
             d(n,2)=d(n,2)-ds2
             d(n,3)=d(n,3)-ds3
             d(n,4)=d(n,4)-ds4
             d(n,5)=d(n,5)-ds5
             d(n-ninc,1)=d(n-ninc,1)+ds1
             d(n-ninc,2)=d(n-ninc,2)+ds2
             d(n-ninc,3)=d(n-ninc,3)+ds3
             d(n-ninc,4)=d(n-ninc,4)+ds4
             d(n-ninc,5)=d(n-ninc,5)+ds5
          enddo
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i1,j1  ,k)
       ind2 = indc(i1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          eps2=ck2*snup(m)
          cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
               sn(m,kdir,2)*sn(m,kdir,2)+ &
               sn(m,kdir,3)*sn(m,kdir,3)
          uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
          vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
          ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
          vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
          rho=v(n,1)+v(n-ninc,1)
          pres=ps(n)+ps(n-ninc)
          a2=gam*pres/rho
!        a2=cson(n)*cson(n-ninc)
          q2=uu**2+vv**2+ww**2
          beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
          rl=0.5*((1.+beta2)*abs(vn) + &
               sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
          get=0.5*(v(n,5)/v(n,1)+v(n-ninc,5)/v(n-ninc,1))
          ge=get-0.5*q2
          gd=(1./beta2-1.)/ge  !matrice de preconditionnement Pc^-1
!
          dw1=v(n,1)-v(n-ninc,1)
          dw2=v(n,2)-v(n-ninc,2)
          dw3=v(n,3)-v(n-ninc,3)
          dw4=v(n,4)-v(n-ninc,4)
          dw5=v(n,5)-v(n-ninc,5)
          prec1=gd*(0.5*q2*dw1-uu*dw2-vv*dw3-ww*dw4+dw5)
          d1= dw1 +    prec1
          d2= dw2 + uu*prec1
          d3= dw3 + vv*prec1
          d4= dw4 + ww*prec1
          d5= dw5 +get*prec1
!
          ds1=eps2*rl*d1
          ds2=eps2*rl*d2
          ds3=eps2*rl*d3
          ds4=eps2*rl*d4
          ds5=eps2*rl*d5
!
          d(n,1)=d(n,1)-ds1
          d(n,2)=d(n,2)-ds2
          d(n,3)=d(n,3)-ds3
          d(n,4)=d(n,4)-ds4
          d(n,5)=d(n,5)-ds5
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i2,j1  ,k)
       ind2 = indc(i2,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          eps2=ck2*snup(m-ninc)
          cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
               sn(m,kdir,2)*sn(m,kdir,2)+ &
               sn(m,kdir,3)*sn(m,kdir,3)
          uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
          vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
          ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
          vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
          rho=v(n,1)+v(n-ninc,1)
          pres=ps(n)+ps(n-ninc)
          a2=gam*pres/rho
!        a2=cson(n)*cson(n-ninc)
          q2=uu**2+vv**2+ww**2
          beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
          rl=0.5*((1.+beta2)*abs(vn) + &
               sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
          get=0.5*(v(n,5)/v(n,1)+v(n-ninc,5)/v(n-ninc,1))
          ge=get-0.5*q2
          gd=(1./beta2-1.)/ge   !matrice de preconditionnement Pc^-1
!
          dw1=v(n,1)-v(n-ninc,1)
          dw2=v(n,2)-v(n-ninc,2)
          dw3=v(n,3)-v(n-ninc,3)
          dw4=v(n,4)-v(n-ninc,4)
          dw5=v(n,5)-v(n-ninc,5)
          prec1=gd*(0.5*q2*dw1-uu*dw2-vv*dw3-ww*dw4+dw5)
          d1= dw1 +    prec1
          d2= dw2 + uu*prec1
          d3= dw3 + vv*prec1
          d4= dw4 + ww*prec1
          d5= dw5 +get*prec1
!
          ds1=eps2*rl*d1
          ds2=eps2*rl*d2
          ds3=eps2*rl*d3
          ds4=eps2*rl*d4
          ds5=eps2*rl*d5
!
          d(n-ninc,1)=d(n-ninc,1)+ds1
          d(n-ninc,2)=d(n-ninc,2)+ds2
          d(n-ninc,3)=d(n-ninc,3)+ds3
          d(n-ninc,4)=d(n-ninc,4)+ds4
          d(n-ninc,5)=d(n-ninc,5)+ds5
       enddo
    enddo
!
!------direction j----------------------------------------------
!
    kdir=2
    ninc=ncj
!
!c---senseur dissipation artificielle
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             snup(m)=abs(ps(n+ninc)-2.*ps(n)+ps(n-ninc)) &
                  /(ps(n+ninc)+2.*ps(n)+ps(n-ninc))
          enddo
       enddo
    enddo
!
    do k=k1,k2m1
       do j=j1p1,j2m1
          ind1 = indc(i1,j,k)
          ind2 = indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             eps2=ck2*max(snup(m-ninc),snup(m))
             eps4=max(0.,ck4-eps2)
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
             a2=cson(n)*cson(n-ninc)
             q2=uu**2+vv**2+ww**2
             beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
             rl=0.5*((1.+beta2)*abs(vn) + &
                  sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
             get=0.5*(v(n,5)/v(n,1)+v(n-ninc,5)/v(n-ninc,1))
             ge=get-0.5*q2
             gd=(1./beta2-1.)/ge    !matrice de preconditionnement Pc^-1
!
             dw1=v(n,1)-v(n-ninc,1)
             dw2=v(n,2)-v(n-ninc,2)
             dw3=v(n,3)-v(n-ninc,3)
             dw4=v(n,4)-v(n-ninc,4)
             dw5=v(n,5)-v(n-ninc,5)
             d3w1=v(n+ninc,1)-3.*v(n,1)+3.*v(n-ninc,1)-v(n-2*ninc,1)
             d3w2=v(n+ninc,2)-3.*v(n,2)+3.*v(n-ninc,2)-v(n-2*ninc,2)
             d3w3=v(n+ninc,3)-3.*v(n,3)+3.*v(n-ninc,3)-v(n-2*ninc,3)
             d3w4=v(n+ninc,4)-3.*v(n,4)+3.*v(n-ninc,4)-v(n-2*ninc,4)
             d3w5=v(n+ninc,5)-3.*v(n,5)+3.*v(n-ninc,5)-v(n-2*ninc,5)
!
             prec1=gd*(0.5*q2*dw1-uu*dw2-vv*dw3-ww*dw4+dw5)
             d1= dw1 +    prec1
             d2= dw2 + uu*prec1
             d3= dw3 + vv*prec1
             d4= dw4 + ww*prec1
             d5= dw5 +get*prec1
             prec3=gd*(0.5*q2*d3w1-uu*d3w2-vv*d3w3-ww*d3w4+d3w5)
             dd1= d3w1 +    prec3
             dd2= d3w2 + uu*prec3
             dd3= d3w3 + vv*prec3
             dd4= d3w4 + ww*prec3
             dd5= d3w5 +get*prec3
!
             ds1=eps2*rl*d1 - eps4*rl*dd1
             ds2=eps2*rl*d2 - eps4*rl*dd2
             ds3=eps2*rl*d3 - eps4*rl*dd3
             ds4=eps2*rl*d4 - eps4*rl*dd4
             ds5=eps2*rl*d5 - eps4*rl*dd5
!
             d(n,1)=d(n,1)-ds1
             d(n,2)=d(n,2)-ds2
             d(n,3)=d(n,3)-ds3
             d(n,4)=d(n,4)-ds4
             d(n,5)=d(n,5)-ds5
             d(n-ninc,1)=d(n-ninc,1)+ds1
             d(n-ninc,2)=d(n-ninc,2)+ds2
             d(n-ninc,3)=d(n-ninc,3)+ds3
             d(n-ninc,4)=d(n-ninc,4)+ds4
             d(n-ninc,5)=d(n-ninc,5)+ds5
          enddo
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j1,k)
       ind2 = indc(i2m1,j1,k)
!$OMP SIMD
       do n=ind1,ind2
          m=n-n0c
          eps2=ck2*snup(m)
          cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
               sn(m,kdir,2)*sn(m,kdir,2)+ &
               sn(m,kdir,3)*sn(m,kdir,3)
          uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
          vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
          ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
          vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
          rho=v(n,1)+v(n-ninc,1)
          pres=ps(n)+ps(n-ninc)
          a2=gam*pres/rho
!         a2=cson(n)*cson(n-ninc)
          q2=uu**2+vv**2+ww**2
          beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
          rl=0.5*((1.+beta2)*abs(vn) + &
               sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
          get=0.5*(v(n,5)/v(n,1)+v(n-ninc,5)/v(n-ninc,1))
          ge=get-0.5*q2
          gd=(1./beta2-1.)/ge    !matrice de preconditionnement Pc^-1
!
          dw1=v(n,1)-v(n-ninc,1)
          dw2=v(n,2)-v(n-ninc,2)
          dw3=v(n,3)-v(n-ninc,3)
          dw4=v(n,4)-v(n-ninc,4)
          dw5=v(n,5)-v(n-ninc,5)
          prec1=gd*(0.5*q2*dw1-uu*dw2-vv*dw3-ww*dw4+dw5)
          d1= dw1 +    prec1
          d2= dw2 + uu*prec1
          d3= dw3 + vv*prec1
          d4= dw4 + ww*prec1
          d5= dw5 +get*prec1
!
          ds1=eps2*rl*d1
          ds2=eps2*rl*d2
          ds3=eps2*rl*d3
          ds4=eps2*rl*d4
          ds5=eps2*rl*d5
!
          d(n,1)=d(n,1)-ds1
          d(n,2)=d(n,2)-ds2
          d(n,3)=d(n,3)-ds3
          d(n,4)=d(n,4)-ds4
          d(n,5)=d(n,5)-ds5
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j2,k)
       ind2 = indc(i2m1,j2,k)
!$OMP SIMD
       do n=ind1,ind2
          m=n-n0c
          eps2=ck2*snup(m-ninc)
          cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
               sn(m,kdir,2)*sn(m,kdir,2)+ &
               sn(m,kdir,3)*sn(m,kdir,3)
          uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
          vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
          ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
          vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
          rho=v(n,1)+v(n-ninc,1)
          pres=ps(n)+ps(n-ninc)
          a2=gam*pres/rho
!         a2=cson(n)*cson(n-ninc)
          q2=uu**2+vv**2+ww**2
          beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
          rl=0.5*((1.+beta2)*abs(vn) + &
               sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
          get=0.5*(v(n,5)/v(n,1)+v(n-ninc,5)/v(n-ninc,1))
          ge=get-0.5*q2
          gd=(1./beta2-1.)/ge    ! matrice de preconditionnement Pc^-1
!
          dw1=v(n,1)-v(n-ninc,1)
          dw2=v(n,2)-v(n-ninc,2)
          dw3=v(n,3)-v(n-ninc,3)
          dw4=v(n,4)-v(n-ninc,4)
          dw5=v(n,5)-v(n-ninc,5)
          prec1=gd*(0.5*q2*dw1-uu*dw2-vv*dw3-ww*dw4+dw5)
          d1= dw1 +    prec1
          d2= dw2 + uu*prec1
          d3= dw3 + vv*prec1
          d4= dw4 + ww*prec1
          d5= dw5 +get*prec1
!
          ds1=eps2*rl*d1
          ds2=eps2*rl*d2
          ds3=eps2*rl*d3
          ds4=eps2*rl*d4
          ds5=eps2*rl*d5
!
          d(n-ninc,1)=d(n-ninc,1)+ds1
          d(n-ninc,2)=d(n-ninc,2)+ds2
          d(n-ninc,3)=d(n-ninc,3)+ds3
          d(n-ninc,4)=d(n-ninc,4)+ds4
          d(n-ninc,5)=d(n-ninc,5)+ds5
       enddo
    enddo
!
!------direction k----------------------------------------------
!
    if(equat(3:4).eq.'3d') then
       kdir=3
       ninc=nck
!
!c---senseur dissipation artificielle--------------------------
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1,j,k)
             ind2 = indc(i2m1,j,k)
!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                snup(m)=abs(ps(n+ninc)-2.*ps(n)+ps(n-ninc)) &
                     /(ps(n+ninc)+2.*ps(n)+ps(n-ninc))
             enddo
          enddo
       enddo
!
       do k=k1p1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                eps2=ck2*max(snup(m-ninc),snup(m))
                eps4=max(0.,ck4-eps2)
                cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                     sn(m,kdir,2)*sn(m,kdir,2)+ &
                     sn(m,kdir,3)*sn(m,kdir,3)
                uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
                vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
                ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
                vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
                a2=cson(n)*cson(n-ninc)
                q2=uu**2+vv**2+ww**2
                beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
                rl=0.5*((1.+beta2)*abs(vn) + &
                     sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
                get=0.5*(v(n,5)/v(n,1)+v(n-ninc,5)/v(n-ninc,1))
                ge=get-0.5*q2
                gd=(1./beta2-1.)/ge   !matrice de preconditionnement Pc^-1
!
                dw1=v(n,1)-v(n-ninc,1)
                dw2=v(n,2)-v(n-ninc,2)
                dw3=v(n,3)-v(n-ninc,3)
                dw4=v(n,4)-v(n-ninc,4)
                dw5=v(n,5)-v(n-ninc,5)
                d3w1=v(n+ninc,1)-3.*v(n,1)+3.*v(n-ninc,1)-v(n-2*ninc,1)
                d3w2=v(n+ninc,2)-3.*v(n,2)+3.*v(n-ninc,2)-v(n-2*ninc,2)
                d3w3=v(n+ninc,3)-3.*v(n,3)+3.*v(n-ninc,3)-v(n-2*ninc,3)
                d3w4=v(n+ninc,4)-3.*v(n,4)+3.*v(n-ninc,4)-v(n-2*ninc,4)
                d3w5=v(n+ninc,5)-3.*v(n,5)+3.*v(n-ninc,5)-v(n-2*ninc,5)
!
                prec1=gd*(0.5*q2*dw1-uu*dw2-vv*dw3-ww*dw4+dw5)
                d1= dw1 +    prec1
                d2= dw2 + uu*prec1
                d3= dw3 + vv*prec1
                d4= dw4 + ww*prec1
                d5= dw5 +get*prec1
                prec3=gd*(0.5*q2*d3w1-uu*d3w2-vv*d3w3-ww*d3w4+d3w5)
                dd1= d3w1 +    prec3
                dd2= d3w2 + uu*prec3
                dd3= d3w3 + vv*prec3
                dd4= d3w4 + ww*prec3
                dd5= d3w5 +get*prec3
!
                ds1=eps2*rl*d1 - eps4*rl*dd1
                ds2=eps2*rl*d2 - eps4*rl*dd2
                ds3=eps2*rl*d3 - eps4*rl*dd3
                ds4=eps2*rl*d4 - eps4*rl*dd4
                ds5=eps2*rl*d5 - eps4*rl*dd5
!
                d(n,1)=d(n,1)-ds1
                d(n,2)=d(n,2)-ds2
                d(n,3)=d(n,3)-ds3
                d(n,4)=d(n,4)-ds4
                d(n,5)=d(n,5)-ds5
                d(n-ninc,1)=d(n-ninc,1)+ds1
                d(n-ninc,2)=d(n-ninc,2)+ds2
                d(n-ninc,3)=d(n-ninc,3)+ds3
                d(n-ninc,4)=d(n-ninc,4)+ds4
                d(n-ninc,5)=d(n-ninc,5)+ds5
             enddo
          enddo
       enddo
!
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k1)
          ind2 = indc(i2m1,j,k1)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             eps2=ck2*snup(m)
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
             rho=v(n,1)+v(n-ninc,1)
             pres=ps(n)+ps(n-ninc)
             a2=gam*pres/rho
!         a2=cson(n)*cson(n-ninc)
             q2=uu**2+vv**2+ww**2
             beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
             rl=0.5*((1.+beta2)*abs(vn) + &
                  sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
             get=0.5*(v(n,5)/v(n,1)+v(n-ninc,5)/v(n-ninc,1))
             ge=get-0.5*q2
             gd=(1./beta2-1.)/ge   !matrice de preconditionnement Pc^-1
!
             dw1=v(n,1)-v(n-ninc,1)
             dw2=v(n,2)-v(n-ninc,2)
             dw3=v(n,3)-v(n-ninc,3)
             dw4=v(n,4)-v(n-ninc,4)
             dw5=v(n,5)-v(n-ninc,5)
             prec1=gd*(0.5*q2*dw1-uu*dw2-vv*dw3-ww*dw4+dw5)
             d1= dw1 +    prec1
             d2= dw2 + uu*prec1
             d3= dw3 + vv*prec1
             d4= dw4 + ww*prec1
             d5= dw5 +get*prec1
!
             ds1=eps2*rl*d1
             ds2=eps2*rl*d2
             ds3=eps2*rl*d3
             ds4=eps2*rl*d4
             ds5=eps2*rl*d5
!
             d(n,1)=d(n,1)-ds1
             d(n,2)=d(n,2)-ds2
             d(n,3)=d(n,3)-ds3
             d(n,4)=d(n,4)-ds4
             d(n,5)=d(n,5)-ds5
          enddo
       enddo
!
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k2)
          ind2 = indc(i2m1,j,k2)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             eps2=ck2*snup(m-ninc)
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
             rho=v(n,1)+v(n-ninc,1)
             pres=ps(n)+ps(n-ninc)
             a2=gam*pres/rho
!         a2=cson(n)*cson(n-ninc)
             q2=uu**2+vv**2+ww**2
             beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
             rl=0.5*((1.+beta2)*abs(vn) + &
                  sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
             get=0.5*(v(n,5)/v(n,1)+v(n-ninc,5)/v(n-ninc,1))
             ge=get-0.5*q2
             gd=(1./beta2-1.)/ge     !matrice de preconditionnement Pc^-1
!
             dw1=v(n,1)-v(n-ninc,1)
             dw2=v(n,2)-v(n-ninc,2)
             dw3=v(n,3)-v(n-ninc,3)
             dw4=v(n,4)-v(n-ninc,4)
             dw5=v(n,5)-v(n-ninc,5)
             prec1=gd*(0.5*q2*dw1-uu*dw2-vv*dw3-ww*dw4+dw5)
             d1= dw1 +    prec1
             d2= dw2 + uu*prec1
             d3= dw3 + vv*prec1
             d4= dw4 + ww*prec1
             d5= dw5 +get*prec1
!
             ds1=eps2*rl*d1
             ds2=eps2*rl*d2
             ds3=eps2*rl*d3
             ds4=eps2*rl*d4
             ds5=eps2*rl*d5
!
             d(n-ninc,1)=d(n-ninc,1)+ds1
             d(n-ninc,2)=d(n-ninc,2)+ds2
             d(n-ninc,3)=d(n-ninc,3)+ds3
             d(n-ninc,4)=d(n-ninc,4)+ds4
             d(n-ninc,5)=d(n-ninc,5)+ds5
          enddo
       enddo
!
    endif
    deallocate(temp)
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
  end subroutine dissip_jameson_prcd2
end module mod_dissip_jameson_prcd2
