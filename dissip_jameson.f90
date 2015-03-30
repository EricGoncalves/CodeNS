module mod_dissip_jameson
  implicit none
contains
  subroutine dissip_jameson( &
       lm,v,d, &
       equat, &
       sn,lgsnlt, &
       snu,ps,cson)
!
!***********************************************************************
!
!_DA  DATE_C : mars 2002 - Eric Goncalves / Sinumef
!
!     ACT
!      Calcul de la dissipation scalaire du schema de Jameson.
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
    integer          ::      i,    i1,  i1m1,  i1p1,    i2
    integer          ::   i2m1,    id,  ind1,  ind2,     j
    integer          ::     j1,  j1m1,  j1p1,    j2,  j2m1
    integer          ::     jd,     k,    k1,  k1m1,  k1p1
    integer          ::     k2,  k2m1,    kd,  kdir,lgsnlt
    integer          ::     lm,     m,     n,   n0c,   nci
    integer          ::    ncj,   nck,   nid,  nijd,  ninc
    integer          ::    njd
    double precision ::                   cc,                 ck2,                 ck4,                cnds,          cson(ip11)
    double precision ::         d(ip11,ip60),                 ds1,                 ds2,                 ds3,                 ds4
    double precision ::                  ds5,                eps2,                eps4,            ps(ip11),                  rl
    double precision :: sn(lgsnlt,nind,ndir),           snu(ip00),                  uu,        v(ip11,ip60),                  vv
    double precision ::                   ww
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
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
!     coefficient de dissipation artificielle champ moyen
    ck2=ki2(lm)
    ck4=ki4(lm)
!
!-----calcul des densites de flux visqueuses----------------------------------------
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
!       ps(n)=gam1*(v(n,5)-pinfl-0.5*(v(n,2)**2+v(n,3)**2+v(n,4)**2)/v(n,1))
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
             snu(m)=abs(ps(n+ninc)-2.*ps(n)+ps(n-ninc)) &
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
             eps2=ck2*max(snu(m-ninc),snu(m))
             eps4=max(0.,ck4-eps2)
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             cc=0.5*(cson(n)+cson(n-ninc))
             rl=abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))+sqrt(cnds)*cc
             ds1=eps2*rl*(v(n,1)-v(n-ninc,1)) &
                  -eps4*rl*(v(n+ninc,1)-3.*v(n,1)+3.*v(n-ninc,1)-v(n-2*ninc,1))
             ds2=eps2*rl*(v(n,2)-v(n-ninc,2)) &
                  -eps4*rl*(v(n+ninc,2)-3.*v(n,2)+3.*v(n-ninc,2)-v(n-2*ninc,2))
             ds3=eps2*rl*(v(n,3)-v(n-ninc,3)) &
                  -eps4*rl*(v(n+ninc,3)-3.*v(n,3)+3.*v(n-ninc,3)-v(n-2*ninc,3))
             ds4=eps2*rl*(v(n,4)-v(n-ninc,4)) &
                  -eps4*rl*(v(n+ninc,4)-3.*v(n,4)+3.*v(n-ninc,4)-v(n-2*ninc,4))
             ds5=eps2*rl*(v(n,5)-v(n-ninc,5)) &
                  -eps4*rl*(v(n+ninc,5)-3.*v(n,5)+3.*v(n-ninc,5)-v(n-2*ninc,5))
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
          eps2=ck2*snu(m)
          cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
               sn(m,kdir,2)*sn(m,kdir,2)+ &
               sn(m,kdir,3)*sn(m,kdir,3)
          uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
          vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
          ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
          cc=0.5*(cson(n)+cson(n-ninc))
          rl=abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))+sqrt(cnds)*cc
          ds1=eps2*rl*(v(n,1)-v(n-ninc,1))
          ds2=eps2*rl*(v(n,2)-v(n-ninc,2))
          ds3=eps2*rl*(v(n,3)-v(n-ninc,3))
          ds4=eps2*rl*(v(n,4)-v(n-ninc,4))
          ds5=eps2*rl*(v(n,5)-v(n-ninc,5))
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
          eps2=ck2*snu(m-ninc)
          cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
               sn(m,kdir,2)*sn(m,kdir,2)+ &
               sn(m,kdir,3)*sn(m,kdir,3)
          uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
          vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
          ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
          cc=0.5*(cson(n)+cson(n-ninc))
          rl=abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))+sqrt(cnds)*cc
          ds1=eps2*rl*(v(n,1)-v(n-ninc,1))
          ds2=eps2*rl*(v(n,2)-v(n-ninc,2))
          ds3=eps2*rl*(v(n,3)-v(n-ninc,3))
          ds4=eps2*rl*(v(n,4)-v(n-ninc,4))
          ds5=eps2*rl*(v(n,5)-v(n-ninc,5))
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
             snu(m)=abs(ps(n+ninc)-2.*ps(n)+ps(n-ninc)) &
                  /(ps(n+ninc)+2.*ps(n)+ps(n-ninc))
          enddo
       enddo
    enddo
!
    do k=k1,k2m1
       do j=j1p1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             eps2=ck2*max(snu(m-ninc),snu(m))
             eps4=max(0.,ck4-eps2)
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             cc=0.5*(cson(n)+cson(n-ninc))
             rl=abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))+sqrt(cnds)*cc
             ds1=eps2*rl*(v(n,1)-v(n-ninc,1)) &
                  -eps4*rl*(v(n+ninc,1)-3.*v(n,1)+3.*v(n-ninc,1)-v(n-2*ninc,1))
             ds2=eps2*rl*(v(n,2)-v(n-ninc,2)) &
                  -eps4*rl*(v(n+ninc,2)-3.*v(n,2)+3.*v(n-ninc,2)-v(n-2*ninc,2))
             ds3=eps2*rl*(v(n,3)-v(n-ninc,3)) &
                  -eps4*rl*(v(n+ninc,3)-3.*v(n,3)+3.*v(n-ninc,3)-v(n-2*ninc,3))
             ds4=eps2*rl*(v(n,4)-v(n-ninc,4)) &
                  -eps4*rl*(v(n+ninc,4)-3.*v(n,4)+3.*v(n-ninc,4)-v(n-2*ninc,4))
             ds5=eps2*rl*(v(n,5)-v(n-ninc,5)) &
                  -eps4*rl*(v(n+ninc,5)-3.*v(n,5)+3.*v(n-ninc,5)-v(n-2*ninc,5))
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
          eps2=ck2*snu(m)
          cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
               sn(m,kdir,2)*sn(m,kdir,2)+ &
               sn(m,kdir,3)*sn(m,kdir,3)
          uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
          vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
          ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
          cc=0.5*(cson(n)+cson(n-ninc))
          rl=abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))+sqrt(cnds)*cc
          ds1=eps2*rl*(v(n,1)-v(n-ninc,1))
          ds2=eps2*rl*(v(n,2)-v(n-ninc,2))
          ds3=eps2*rl*(v(n,3)-v(n-ninc,3))
          ds4=eps2*rl*(v(n,4)-v(n-ninc,4))
          ds5=eps2*rl*(v(n,5)-v(n-ninc,5))
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
          eps2=ck2*snu(m-ninc)
          cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
               sn(m,kdir,2)*sn(m,kdir,2)+ &
               sn(m,kdir,3)*sn(m,kdir,3)
          uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
          vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
          ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
          cc=0.5*(cson(n)+cson(n-ninc))
          rl=abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))+sqrt(cnds)*cc
          ds1=eps2*rl*(v(n,1)-v(n-ninc,1))
          ds2=eps2*rl*(v(n,2)-v(n-ninc,2))
          ds3=eps2*rl*(v(n,3)-v(n-ninc,3))
          ds4=eps2*rl*(v(n,4)-v(n-ninc,4))
          ds5=eps2*rl*(v(n,5)-v(n-ninc,5))
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
!----senseur dissipation artificielle
       do k=k1,k2m1
          do j=j1,j2m1
             ind1=indc(i1  ,j,k)
             ind2=indc(i2m1,j,k)
!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                snu(m)=abs(ps(n+ninc)-2.*ps(n)+ps(n-ninc)) &
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
                eps2=ck2*max(snu(m-ninc),snu(m))
                eps4=max(0.,ck4-eps2)
                cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                     sn(m,kdir,2)*sn(m,kdir,2)+ &
                     sn(m,kdir,3)*sn(m,kdir,3)
                uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
                vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
                ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
                cc=0.5*(cson(n)+cson(n-ninc))
                rl=abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))+sqrt(cnds)*cc
                ds1=eps2*rl*(v(n,1)-v(n-ninc,1)) &
                     -eps4*rl*(v(n+ninc,1)-3.*v(n,1)+3.*v(n-ninc,1)-v(n-2*ninc,1))
                ds2=eps2*rl*(v(n,2)-v(n-ninc,2)) &
                     -eps4*rl*(v(n+ninc,2)-3.*v(n,2)+3.*v(n-ninc,2)-v(n-2*ninc,2))
                ds3=eps2*rl*(v(n,3)-v(n-ninc,3)) &
                     -eps4*rl*(v(n+ninc,3)-3.*v(n,3)+3.*v(n-ninc,3)-v(n-2*ninc,3))
                ds4=eps2*rl*(v(n,4)-v(n-ninc,4)) &
                     -eps4*rl*(v(n+ninc,4)-3.*v(n,4)+3.*v(n-ninc,4)-v(n-2*ninc,4))
                ds5=eps2*rl*(v(n,5)-v(n-ninc,5)) &
                     -eps4*rl*(v(n+ninc,5)-3.*v(n,5)+3.*v(n-ninc,5)-v(n-2*ninc,5))
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
             eps2=ck2*snu(m)
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             cc=0.5*(cson(n)+cson(n-ninc))
             rl=abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))+sqrt(cnds)*cc
             ds1=eps2*rl*(v(n,1)-v(n-ninc,1))
             ds2=eps2*rl*(v(n,2)-v(n-ninc,2))
             ds3=eps2*rl*(v(n,3)-v(n-ninc,3))
             ds4=eps2*rl*(v(n,4)-v(n-ninc,4))
             ds5=eps2*rl*(v(n,5)-v(n-ninc,5))
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
             eps2=ck2*snu(m-ninc)
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             cc=0.5*(cson(n)+cson(n-ninc))
             rl=abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))+sqrt(cnds)*cc
             ds1=eps2*rl*(v(n,1)-v(n-ninc,1))
             ds2=eps2*rl*(v(n,2)-v(n-ninc,2))
             ds3=eps2*rl*(v(n,3)-v(n-ninc,3))
             ds4=eps2*rl*(v(n,4)-v(n-ninc,4))
             ds5=eps2*rl*(v(n,5)-v(n-ninc,5))
             d(n-ninc,1)=d(n-ninc,1)+ds1
             d(n-ninc,2)=d(n-ninc,2)+ds2
             d(n-ninc,3)=d(n-ninc,3)+ds3
             d(n-ninc,4)=d(n-ninc,4)+ds4
             d(n-ninc,5)=d(n-ninc,5)+ds5
          enddo
       enddo
!
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
  end subroutine dissip_jameson
end module mod_dissip_jameson
