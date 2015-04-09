module mod_sch_jameson3
  implicit none
contains
  subroutine sch_jameson3( &
       lm,ityprk, &
       u,v,d,ff, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       equat, &
       sn,lgsnlt, &
       fxx,fyy,fzz,fxy,fxz,fyz,fex,fey,fez, &
       ps)
!
!***********************************************************************
!
!_DA  DATE_C : avril 2002 - Eric Goncalves / Sinumef
!
!     ACT
!_A    Calcul des bilans de flux physiques pour chaque cellule.
!_A    Schema de Jameson avec correction de l'erreur dispersive.
!_A    Ordre 3 en cartesien.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    implicit none
    integer          ::      i,    i1,  i1m1,  i1p1,    i2
    integer          ::   i2m1,    id,  ind1,  ind2,ityprk
    integer          ::      j,    j1,  j1m1,  j1p1,    j2
    integer          ::   j2m1,    jd,     k,    k1,  k1m1
    integer          ::   k1p1,    k2,  k2m1,    kd,  kdir
    integer          :: lgsnlt,    lm,     m,    m2,     n
    integer          ::    n0c,    n2,   nci,   ncj,   nck
    integer          ::    nid,  nijd,  ninc,   njd
    double precision ::                   c0,        d(ip11,ip60),           fex(ip00),           fey(ip00),           fez(ip00)
    double precision ::        ff(ip11,ip60),           fxx(ip00),           fxy(ip00),           fxz(ip00),           fyy(ip00)
    double precision ::            fyz(ip00),           fzz(ip00),            ps(ip11),           qcx(ip12),           qcy(ip12)
    double precision ::            qcz(ip12),                 si0,                 si1,                 si2,                 si3
    double precision ::                  si4,                 sj0,                 sj1,                 sj2,                 sj3
    double precision ::                  sj4,                 sk0,                 sk1,                 sk2,                 sk3
    double precision ::                  sk4,sn(lgsnlt,nind,ndir),          toxx(ip12),          toxy(ip12),          toxz(ip12)
    double precision ::           toyy(ip12),          toyz(ip12),          tozz(ip12),        u(ip11,ip60),        v(ip11,ip60),norm
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
!-----calcul des densites de flux visqueuses--------------------------------
!
    if(equat(3:5).eq.'2dk') then
       ind1 = indc(i1m1,j1m1,k1  )
       ind2 = indc(i2  ,j2  ,k2m1)
    elseif(equat(3:4).eq.'3d') then
       ind1 = indc(i1m1,j1m1,k1m1)
       ind2 = indc(i2  ,j2  ,k2  )
    endif
!
!!$OMP DO SIMD
!$OMP SIMD
       do n=ind1,ind2
          m=n-n0c
          u(n,1)=0.
          u(n,2)=0.
          u(n,3)=0.
          u(n,4)=0.
          u(n,5)=0.
        norm=1./v(n,1)
        fxx(m)=v(n,2)*(v(n,2)*norm)+ps(n)-pinfl
        fxy(m)=v(n,3)*(v(n,2)*norm)
        fxz(m)=v(n,4)*(v(n,2)*norm)
        fyy(m)=v(n,3)*(v(n,3)*norm)+ps(n)-pinfl
        fyz(m)=v(n,4)*(v(n,3)*norm)
        fzz(m)=v(n,4)*(v(n,4)*norm)+ps(n)-pinfl
        fex(m)=(v(n,5)+ps(n)-pinfl)*v(n,2)*norm
        fey(m)=(v(n,5)+ps(n)-pinfl)*v(n,3)*norm
        fez(m)=(v(n,5)+ps(n)-pinfl)*v(n,4)*norm
       enddo
!!$OMP END DO SIMD
    if (equat(1:2).eq.'ns') then
!!$OMP DO SIMD
!$OMP SIMD
       do n=ind1,ind2
          m=n-n0c
          fxx(m)=fxx(m)-toxx(n)
          fxy(m)=fxy(m)-toxy(n)
          fxz(m)=fxz(m)-toxz(n)
          fyy(m)=fyy(m)-toyy(n)
          fyz(m)=fyz(m)-toyz(n)
          fzz(m)=fzz(m)-tozz(n)
          fex(m)=fex(m)-(toxx(n)*v(n,2)+toxy(n)*v(n,3)+toxz(n)*v(n,4))/v(n,1)-qcx(n)
          fey(m)=fey(m)-(toyy(n)*v(n,3)+toxy(n)*v(n,2)+toyz(n)*v(n,4))/v(n,1)-qcy(n)
          fez(m)=fez(m)-(tozz(n)*v(n,4)+toxz(n)*v(n,2)+toyz(n)*v(n,3))/v(n,1)-qcz(n)
       enddo
!!$OMP END DO SIMD
    endif
!
!     coefficient
    c0=1./6.
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
!!$OMP DO collapse(2)
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1+2,j,k)
          ind2 = indc(i2-2,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             n2=n-2*ninc
             m2=m-2*ninc
             si0=c0*((7.*(v(n,2)+v(n-ninc,2))-v(n2,2)-v(n+ninc,2))*sn(m,kdir,1) &
                  +(7.*(v(n,3)+v(n-ninc,3))-v(n2,3)-v(n+ninc,3))*sn(m,kdir,2) &
                  +(7.*(v(n,4)+v(n-ninc,4))-v(n2,4)-v(n+ninc,4))*sn(m,kdir,3))
             si1=c0*((7.*(fxx(m)+fxx(m-ninc))-fxx(m2)-fxx(m+ninc))*sn(m,kdir,1) &
                  +(7.*(fxy(m)+fxy(m-ninc))-fxy(m2)-fxy(m+ninc))*sn(m,kdir,2) &
                  +(7.*(fxz(m)+fxz(m-ninc))-fxz(m2)-fxz(m+ninc))*sn(m,kdir,3))
             si2=c0*((7.*(fxy(m)+fxy(m-ninc))-fxy(m2)-fxy(m+ninc))*sn(m,kdir,1) &
                  +(7.*(fyy(m)+fyy(m-ninc))-fyy(m2)-fyy(m+ninc))*sn(m,kdir,2) &
                  +(7.*(fyz(m)+fyz(m-ninc))-fyz(m2)-fyz(m+ninc))*sn(m,kdir,3))
             si3=c0*((7.*(fxz(m)+fxz(m-ninc))-fxz(m2)-fxz(m+ninc))*sn(m,kdir,1) &
                  +(7.*(fyz(m)+fyz(m-ninc))-fyz(m2)-fyz(m+ninc))*sn(m,kdir,2) &
                  +(7.*(fzz(m)+fzz(m-ninc))-fzz(m2)-fzz(m+ninc))*sn(m,kdir,3))
             si4=c0*((7.*(fex(m)+fex(m-ninc))-fex(m2)-fex(m+ninc))*sn(m,kdir,1) &
                  +(7.*(fey(m)+fey(m-ninc))-fey(m2)-fey(m+ninc))*sn(m,kdir,2) &
                  +(7.*(fez(m)+fez(m-ninc))-fez(m2)-fez(m+ninc))*sn(m,kdir,3))
             u(n,1)=u(n,1)-si0
             u(n,2)=u(n,2)-si1
             u(n,3)=u(n,3)-si2
             u(n,4)=u(n,4)-si3
             u(n,5)=u(n,5)-si4
             u(n-ninc,1)=u(n-ninc,1)+si0
             u(n-ninc,2)=u(n-ninc,2)+si1
             u(n-ninc,3)=u(n-ninc,3)+si2
             u(n-ninc,4)=u(n-ninc,4)+si3
             u(n-ninc,5)=u(n-ninc,5)+si4
          enddo
       enddo
    enddo
!!$OMP END DO
!
!!$OMP DO
    do k=k1,k2m1
       ind1 = indc(i1p1,j1  ,k)
       ind2 = indc(i1p1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          si0= (v(n,2)+v(n-ninc,2))*sn(m,kdir,1) &
               +(v(n,3)+v(n-ninc,3))*sn(m,kdir,2) &
               +(v(n,4)+v(n-ninc,4))*sn(m,kdir,3)
          si1= (fxx(m)+fxx(m-ninc))*sn(m,kdir,1) &
               +(fxy(m)+fxy(m-ninc))*sn(m,kdir,2) &
               +(fxz(m)+fxz(m-ninc))*sn(m,kdir,3)
          si2= (fxy(m)+fxy(m-ninc))*sn(m,kdir,1) &
               +(fyy(m)+fyy(m-ninc))*sn(m,kdir,2) &
               +(fyz(m)+fyz(m-ninc))*sn(m,kdir,3)
          si3= (fxz(m)+fxz(m-ninc))*sn(m,kdir,1) &
               +(fyz(m)+fyz(m-ninc))*sn(m,kdir,2) &
               +(fzz(m)+fzz(m-ninc))*sn(m,kdir,3)
          si4= (fex(m)+fex(m-ninc))*sn(m,kdir,1) &
               +(fey(m)+fey(m-ninc))*sn(m,kdir,2) &
               +(fez(m)+fez(m-ninc))*sn(m,kdir,3)
          u(n,1)=u(n,1)-si0
          u(n,2)=u(n,2)-si1
          u(n,3)=u(n,3)-si2
          u(n,4)=u(n,4)-si3
          u(n,5)=u(n,5)-si4
          u(n-ninc,1)=u(n-ninc,1)+si0
          u(n-ninc,2)=u(n-ninc,2)+si1
          u(n-ninc,3)=u(n-ninc,3)+si2
          u(n-ninc,4)=u(n-ninc,4)+si3
          u(n-ninc,5)=u(n-ninc,5)+si4
       enddo
    enddo
!!$OMP END DO
!
!!$OMP DO
    do k=k1,k2m1
       ind1 = indc(i2m1,j1  ,k)
       ind2 = indc(i2m1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          si0= (v(n,2)+v(n-ninc,2))*sn(m,kdir,1) &
               +(v(n,3)+v(n-ninc,3))*sn(m,kdir,2) &
               +(v(n,4)+v(n-ninc,4))*sn(m,kdir,3)
          si1= (fxx(m)+fxx(m-ninc))*sn(m,kdir,1) &
               +(fxy(m)+fxy(m-ninc))*sn(m,kdir,2) &
               +(fxz(m)+fxz(m-ninc))*sn(m,kdir,3)
          si2= (fxy(m)+fxy(m-ninc))*sn(m,kdir,1) &
               +(fyy(m)+fyy(m-ninc))*sn(m,kdir,2) &
               +(fyz(m)+fyz(m-ninc))*sn(m,kdir,3)
          si3= (fxz(m)+fxz(m-ninc))*sn(m,kdir,1) &
               +(fyz(m)+fyz(m-ninc))*sn(m,kdir,2) &
               +(fzz(m)+fzz(m-ninc))*sn(m,kdir,3)
          si4= (fex(m)+fex(m-ninc))*sn(m,kdir,1) &
               +(fey(m)+fey(m-ninc))*sn(m,kdir,2) &
               +(fez(m)+fez(m-ninc))*sn(m,kdir,3)
          u(n,1)=u(n,1)-si0
          u(n,2)=u(n,2)-si1
          u(n,3)=u(n,3)-si2
          u(n,4)=u(n,4)-si3
          u(n,5)=u(n,5)-si4
          u(n-ninc,1)=u(n-ninc,1)+si0
          u(n-ninc,2)=u(n-ninc,2)+si1
          u(n-ninc,3)=u(n-ninc,3)+si2
          u(n-ninc,4)=u(n-ninc,4)+si3
          u(n-ninc,5)=u(n-ninc,5)+si4
       enddo
    enddo
!!$OMP END DO
!
!!$OMP DO
    do k=k1,k2m1
       ind1 = indc(i1,j1  ,k)
       ind2 = indc(i1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          si0= 2.*(v(n-ninc,2)*sn(m,kdir,1) &
               +v(n-ninc,3)*sn(m,kdir,2) &
               +v(n-ninc,4)*sn(m,kdir,3))
          si1= 2.*(fxx(m-ninc)*sn(m,kdir,1) &
               +fxy(m-ninc)*sn(m,kdir,2) &
               +fxz(m-ninc)*sn(m,kdir,3))
          si2= 2.*(fxy(m-ninc)*sn(m,kdir,1) &
               +fyy(m-ninc)*sn(m,kdir,2) &
               +fyz(m-ninc)*sn(m,kdir,3))
          si3= 2.*(fxz(m-ninc)*sn(m,kdir,1) &
               +fyz(m-ninc)*sn(m,kdir,2) &
               +fzz(m-ninc)*sn(m,kdir,3))
          si4= 2.*(fex(m-ninc)*sn(m,kdir,1) &
               +fey(m-ninc)*sn(m,kdir,2) &
               +fez(m-ninc)*sn(m,kdir,3))
          u(n,1)=u(n,1)-si0
          u(n,2)=u(n,2)-si1
          u(n,3)=u(n,3)-si2
          u(n,4)=u(n,4)-si3
          u(n,5)=u(n,5)-si4
       enddo
    enddo
!!$OMP END DO
!
!!$OMP DO
    do k=k1,k2m1
       ind1 = indc(i2,j1  ,k)
       ind2 = indc(i2,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          si0= 2.*(v(n,2)*sn(m,kdir,1) &
               +v(n,3)*sn(m,kdir,2) &
               +v(n,4)*sn(m,kdir,3))
          si1= 2.*(fxx(m)*sn(m,kdir,1) &
               +fxy(m)*sn(m,kdir,2) &
               +fxz(m)*sn(m,kdir,3))
          si2= 2.*(fxy(m)*sn(m,kdir,1) &
               +fyy(m)*sn(m,kdir,2) &
               +fyz(m)*sn(m,kdir,3))
          si3= 2.*(fxz(m)*sn(m,kdir,1) &
               +fyz(m)*sn(m,kdir,2) &
               +fzz(m)*sn(m,kdir,3))
          si4= 2.*(fex(m)*sn(m,kdir,1) &
               +fey(m)*sn(m,kdir,2) &
               +fez(m)*sn(m,kdir,3))
          u(n-ninc,1)=u(n-ninc,1)+si0
          u(n-ninc,2)=u(n-ninc,2)+si1
          u(n-ninc,3)=u(n-ninc,3)+si2
          u(n-ninc,4)=u(n-ninc,4)+si3
          u(n-ninc,5)=u(n-ninc,5)+si4
       enddo
    enddo
!!$OMP END DO
!
!------direction j------------------------------------------------------------
!
    kdir=2
    ninc=ncj
!
!!$OMP DO collapse(2)
    do k=k1,k2m1
       do j=j1+2,j2-2
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             n2=n-2*ninc
             m2=m-2*ninc
             sj0=c0*((7.*(v(n,2)+v(n-ninc,2))-v(n2,2)-v(n+ninc,2))*sn(m,kdir,1) &
                  +(7.*(v(n,3)+v(n-ninc,3))-v(n2,3)-v(n+ninc,3))*sn(m,kdir,2) &
                  +(7.*(v(n,4)+v(n-ninc,4))-v(n2,4)-v(n+ninc,4))*sn(m,kdir,3))
             sj1=c0*((7.*(fxx(m)+fxx(m-ninc))-fxx(m2)-fxx(m+ninc))*sn(m,kdir,1) &
                  +(7.*(fxy(m)+fxy(m-ninc))-fxy(m2)-fxy(m+ninc))*sn(m,kdir,2) &
                  +(7.*(fxz(m)+fxz(m-ninc))-fxz(m2)-fxz(m+ninc))*sn(m,kdir,3))
             sj2=c0*((7.*(fxy(m)+fxy(m-ninc))-fxy(m2)-fxy(m+ninc))*sn(m,kdir,1) &
                  +(7.*(fyy(m)+fyy(m-ninc))-fyy(m2)-fyy(m+ninc))*sn(m,kdir,2) &
                  +(7.*(fyz(m)+fyz(m-ninc))-fyz(m2)-fyz(m+ninc))*sn(m,kdir,3))
             sj3=c0*((7.*(fxz(m)+fxz(m-ninc))-fxz(m2)-fxz(m+ninc))*sn(m,kdir,1) &
                  +(7.*(fyz(m)+fyz(m-ninc))-fyz(m2)-fyz(m+ninc))*sn(m,kdir,2) &
                  +(7.*(fzz(m)+fzz(m-ninc))-fzz(m2)-fzz(m+ninc))*sn(m,kdir,3))
             sj4=c0*((7.*(fex(m)+fex(m-ninc))-fex(m2)-fex(m+ninc))*sn(m,kdir,1) &
                  +(7.*(fey(m)+fey(m-ninc))-fey(m2)-fey(m+ninc))*sn(m,kdir,2) &
                  +(7.*(fez(m)+fez(m-ninc))-fez(m2)-fez(m+ninc))*sn(m,kdir,3))
             u(n,1)=u(n,1)-sj0
             u(n,2)=u(n,2)-sj1
             u(n,3)=u(n,3)-sj2
             u(n,4)=u(n,4)-sj3
             u(n,5)=u(n,5)-sj4
             u(n-ninc,1)=u(n-ninc,1)+sj0
             u(n-ninc,2)=u(n-ninc,2)+sj1
             u(n-ninc,3)=u(n-ninc,3)+sj2
             u(n-ninc,4)=u(n-ninc,4)+sj3
             u(n-ninc,5)=u(n-ninc,5)+sj4
          enddo
       enddo
    enddo
!!$OMP END DO
!
!!$OMP DO
    do k=k1,k2m1
       ind1 = indc(i1  ,j1p1,k)
       ind2 = indc(i2m1,j1p1,k)
!$OMP SIMD
       do n=ind1,ind2
          m=n-n0c
          sj0= (v(n,2)+v(n-ninc,2))*sn(m,kdir,1) &
               +(v(n,3)+v(n-ninc,3))*sn(m,kdir,2) &
               +(v(n,4)+v(n-ninc,4))*sn(m,kdir,3)
          sj1= (fxx(m)+fxx(m-ninc))*sn(m,kdir,1) &
               +(fxy(m)+fxy(m-ninc))*sn(m,kdir,2) &
               +(fxz(m)+fxz(m-ninc))*sn(m,kdir,3)
          sj2= (fxy(m)+fxy(m-ninc))*sn(m,kdir,1) &
               +(fyy(m)+fyy(m-ninc))*sn(m,kdir,2) &
               +(fyz(m)+fyz(m-ninc))*sn(m,kdir,3)
          sj3= (fxz(m)+fxz(m-ninc))*sn(m,kdir,1) &
               +(fyz(m)+fyz(m-ninc))*sn(m,kdir,2) &
               +(fzz(m)+fzz(m-ninc))*sn(m,kdir,3)
          sj4= (fex(m)+fex(m-ninc))*sn(m,kdir,1) &
               +(fey(m)+fey(m-ninc))*sn(m,kdir,2) &
               +(fez(m)+fez(m-ninc))*sn(m,kdir,3)
          u(n,1)=u(n,1)-sj0
          u(n,2)=u(n,2)-sj1
          u(n,3)=u(n,3)-sj2
          u(n,4)=u(n,4)-sj3
          u(n,5)=u(n,5)-sj4
          u(n-ninc,1)=u(n-ninc,1)+sj0
          u(n-ninc,2)=u(n-ninc,2)+sj1
          u(n-ninc,3)=u(n-ninc,3)+sj2
          u(n-ninc,4)=u(n-ninc,4)+sj3
          u(n-ninc,5)=u(n-ninc,5)+sj4
       enddo
    enddo
!!$OMP END DO
!
!!$OMP DO
    do k=k1,k2m1
       ind1 = indc(i1  ,j2m1,k)
       ind2 = indc(i2m1,j2m1,k)
!$OMP SIMD
       do n=ind1,ind2
          m=n-n0c
          sj0= (v(n,2)+v(n-ninc,2))*sn(m,kdir,1) &
               +(v(n,3)+v(n-ninc,3))*sn(m,kdir,2) &
               +(v(n,4)+v(n-ninc,4))*sn(m,kdir,3)
          sj1= (fxx(m)+fxx(m-ninc))*sn(m,kdir,1) &
               +(fxy(m)+fxy(m-ninc))*sn(m,kdir,2) &
               +(fxz(m)+fxz(m-ninc))*sn(m,kdir,3)
          sj2= (fxy(m)+fxy(m-ninc))*sn(m,kdir,1) &
               +(fyy(m)+fyy(m-ninc))*sn(m,kdir,2) &
               +(fyz(m)+fyz(m-ninc))*sn(m,kdir,3)
          sj3= (fxz(m)+fxz(m-ninc))*sn(m,kdir,1) &
               +(fyz(m)+fyz(m-ninc))*sn(m,kdir,2) &
               +(fzz(m)+fzz(m-ninc))*sn(m,kdir,3)
          sj4= (fex(m)+fex(m-ninc))*sn(m,kdir,1) &
               +(fey(m)+fey(m-ninc))*sn(m,kdir,2) &
               +(fez(m)+fez(m-ninc))*sn(m,kdir,3)
          u(n,1)=u(n,1)-sj0
          u(n,2)=u(n,2)-sj1
          u(n,3)=u(n,3)-sj2
          u(n,4)=u(n,4)-sj3
          u(n,5)=u(n,5)-sj4
          u(n-ninc,1)=u(n-ninc,1)+sj0
          u(n-ninc,2)=u(n-ninc,2)+sj1
          u(n-ninc,3)=u(n-ninc,3)+sj2
          u(n-ninc,4)=u(n-ninc,4)+sj3
          u(n-ninc,5)=u(n-ninc,5)+sj4
       enddo
    enddo
!!$OMP END DO
!
!!$OMP DO
    do k=k1,k2m1
       ind1 = indc(i1  ,j1,k)
       ind2 = indc(i2m1,j1,k)
!$OMP SIMD
       do n=ind1,ind2
          m=n-n0c
          sj0= 2.*(v(n-ninc,2)*sn(m,kdir,1) &
               +v(n-ninc,3)*sn(m,kdir,2) &
               +v(n-ninc,4)*sn(m,kdir,3))
          sj1= 2.*(fxx(m-ninc)*sn(m,kdir,1) &
               +fxy(m-ninc)*sn(m,kdir,2) &
               +fxz(m-ninc)*sn(m,kdir,3))
          sj2= 2.*(fxy(m-ninc)*sn(m,kdir,1) &
               +fyy(m-ninc)*sn(m,kdir,2) &
               +fyz(m-ninc)*sn(m,kdir,3))
          sj3= 2.*(fxz(m-ninc)*sn(m,kdir,1) &
               +fyz(m-ninc)*sn(m,kdir,2) &
               +fzz(m-ninc)*sn(m,kdir,3))
          sj4= 2.*(fex(m-ninc)*sn(m,kdir,1) &
               +fey(m-ninc)*sn(m,kdir,2) &
               +fez(m-ninc)*sn(m,kdir,3))
          u(n,1)=u(n,1)-sj0
          u(n,2)=u(n,2)-sj1
          u(n,3)=u(n,3)-sj2
          u(n,4)=u(n,4)-sj3
          u(n,5)=u(n,5)-sj4
       enddo
    enddo
!!$OMP END DO
!
!!$OMP DO
    do k=k1,k2m1
       ind1 = indc(i1  ,j2,k)
       ind2 = indc(i2m1,j2,k)
!$OMP SIMD
       do n=ind1,ind2
          m=n-n0c
          sj0= 2.*(v(n,2)*sn(m,kdir,1) &
               +v(n,3)*sn(m,kdir,2) &
               +v(n,4)*sn(m,kdir,3))
          sj1= 2.*(fxx(m)*sn(m,kdir,1) &
               +fxy(m)*sn(m,kdir,2) &
               +fxz(m)*sn(m,kdir,3))
          sj2= 2.*(fxy(m)*sn(m,kdir,1) &
               +fyy(m)*sn(m,kdir,2) &
               +fyz(m)*sn(m,kdir,3))
          sj3= 2.*(fxz(m)*sn(m,kdir,1) &
               +fyz(m)*sn(m,kdir,2) &
               +fzz(m)*sn(m,kdir,3))
          sj4= 2.*(fex(m)*sn(m,kdir,1) &
               +fey(m)*sn(m,kdir,2) &
               +fez(m)*sn(m,kdir,3))
          u(n-ninc,1)=u(n-ninc,1)+sj0
          u(n-ninc,2)=u(n-ninc,2)+sj1
          u(n-ninc,3)=u(n-ninc,3)+sj2
          u(n-ninc,4)=u(n-ninc,4)+sj3
          u(n-ninc,5)=u(n-ninc,5)+sj4
       enddo
    enddo
!!$OMP END DO
!
!------direction k----------------------------------------------
!
    if(equat(3:4).eq.'3d') then
       kdir=3
       ninc=nck
!
!!$OMP DO collapse(2)
       do k=k1+2,k2-2
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                n2=n-2*ninc
                m2=m-2*ninc
                sk0=c0*((7.*(v(n,2)+v(n-ninc,2))-v(n2,2)-v(n+ninc,2))*sn(m,kdir,1) &
                     +(7.*(v(n,3)+v(n-ninc,3))-v(n2,3)-v(n+ninc,3))*sn(m,kdir,2) &
                     +(7.*(v(n,4)+v(n-ninc,4))-v(n2,4)-v(n+ninc,4))*sn(m,kdir,3))
                sk1=c0*((7.*(fxx(m)+fxx(m-ninc))-fxx(m2)-fxx(m+ninc))*sn(m,kdir,1) &
                     +(7.*(fxy(m)+fxy(m-ninc))-fxy(m2)-fxy(m+ninc))*sn(m,kdir,2) &
                     +(7.*(fxz(m)+fxz(m-ninc))-fxz(m2)-fxz(m+ninc))*sn(m,kdir,3))
                sk2=c0*((7.*(fxy(m)+fxy(m-ninc))-fxy(m2)-fxy(m+ninc))*sn(m,kdir,1) &
                     +(7.*(fyy(m)+fyy(m-ninc))-fyy(m2)-fyy(m+ninc))*sn(m,kdir,2) &
                     +(7.*(fyz(m)+fyz(m-ninc))-fyz(m2)-fyz(m+ninc))*sn(m,kdir,3))
                sk3=c0*((7.*(fxz(m)+fxz(m-ninc))-fxz(m2)-fxz(m+ninc))*sn(m,kdir,1) &
                     +(7.*(fyz(m)+fyz(m-ninc))-fyz(m2)-fyz(m+ninc))*sn(m,kdir,2) &
                     +(7.*(fzz(m)+fzz(m-ninc))-fzz(m2)-fzz(m+ninc))*sn(m,kdir,3))
                sk4=c0*((7.*(fex(m)+fex(m-ninc))-fex(m2)-fex(m+ninc))*sn(m,kdir,1) &
                     +(7.*(fey(m)+fey(m-ninc))-fey(m2)-fey(m+ninc))*sn(m,kdir,2) &
                     +(7.*(fez(m)+fez(m-ninc))-fez(m2)-fez(m+ninc))*sn(m,kdir,3))
                u(n,1)=u(n,1)-sk0
                u(n,2)=u(n,2)-sk1
                u(n,3)=u(n,3)-sk2
                u(n,4)=u(n,4)-sk3
                u(n,5)=u(n,5)-sk4
                u(n-ninc,1)=u(n-ninc,1)+sk0
                u(n-ninc,2)=u(n-ninc,2)+sk1
                u(n-ninc,3)=u(n-ninc,3)+sk2
                u(n-ninc,4)=u(n-ninc,4)+sk3
                u(n-ninc,5)=u(n-ninc,5)+sk4
             enddo
          enddo
       enddo
!!$OMP END DO
!
!!$OMP DO
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k1p1)
          ind2 = indc(i2m1,j,k1p1)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             sk0= (v(n,2)+v(n-ninc,2))*sn(m,kdir,1) &
                  +(v(n,3)+v(n-ninc,3))*sn(m,kdir,2) &
                  +(v(n,4)+v(n-ninc,4))*sn(m,kdir,3)
             sk1= (fxx(m)+fxx(m-ninc))*sn(m,kdir,1) &
                  +(fxy(m)+fxy(m-ninc))*sn(m,kdir,2) &
                  +(fxz(m)+fxz(m-ninc))*sn(m,kdir,3)
             sk2= (fxy(m)+fxy(m-ninc))*sn(m,kdir,1) &
                  +(fyy(m)+fyy(m-ninc))*sn(m,kdir,2) &
                  +(fyz(m)+fyz(m-ninc))*sn(m,kdir,3)
             sk3= (fxz(m)+fxz(m-ninc))*sn(m,kdir,1) &
                  +(fyz(m)+fyz(m-ninc))*sn(m,kdir,2) &
                  +(fzz(m)+fzz(m-ninc))*sn(m,kdir,3)
             sk4= (fex(m)+fex(m-ninc))*sn(m,kdir,1) &
                  +(fey(m)+fey(m-ninc))*sn(m,kdir,2) &
                  +(fez(m)+fez(m-ninc))*sn(m,kdir,3)
             u(n,1)=u(n,1)-sk0
             u(n,2)=u(n,2)-sk1
             u(n,3)=u(n,3)-sk2
             u(n,4)=u(n,4)-sk3
             u(n,5)=u(n,5)-sk4
             u(n-ninc,1)=u(n-ninc,1)+sk0
             u(n-ninc,2)=u(n-ninc,2)+sk1
             u(n-ninc,3)=u(n-ninc,3)+sk2
             u(n-ninc,4)=u(n-ninc,4)+sk3
             u(n-ninc,5)=u(n-ninc,5)+sk4
          enddo
       enddo
!!$OMP END DO
!
!!$OMP DO
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k2m1)
          ind2 = indc(i2m1,j,k2m1)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             sk0= (v(n,2)+v(n-ninc,2))*sn(m,kdir,1) &
                  +(v(n,3)+v(n-ninc,3))*sn(m,kdir,2) &
                  +(v(n,4)+v(n-ninc,4))*sn(m,kdir,3)
             sk1= (fxx(m)+fxx(m-ninc))*sn(m,kdir,1) &
                  +(fxy(m)+fxy(m-ninc))*sn(m,kdir,2) &
                  +(fxz(m)+fxz(m-ninc))*sn(m,kdir,3)
             sk2= (fxy(m)+fxy(m-ninc))*sn(m,kdir,1) &
                  +(fyy(m)+fyy(m-ninc))*sn(m,kdir,2) &
                  +(fyz(m)+fyz(m-ninc))*sn(m,kdir,3)
             sk3= (fxz(m)+fxz(m-ninc))*sn(m,kdir,1) &
                  +(fyz(m)+fyz(m-ninc))*sn(m,kdir,2) &
                  +(fzz(m)+fzz(m-ninc))*sn(m,kdir,3)
             sk4= (fex(m)+fex(m-ninc))*sn(m,kdir,1) &
                  +(fey(m)+fey(m-ninc))*sn(m,kdir,2) &
                  +(fez(m)+fez(m-ninc))*sn(m,kdir,3)
             u(n,1)=u(n,1)-sk0
             u(n,2)=u(n,2)-sk1
             u(n,3)=u(n,3)-sk2
             u(n,4)=u(n,4)-sk3
             u(n,5)=u(n,5)-sk4
             u(n-ninc,1)=u(n-ninc,1)+sk0
             u(n-ninc,2)=u(n-ninc,2)+sk1
             u(n-ninc,3)=u(n-ninc,3)+sk2
             u(n-ninc,4)=u(n-ninc,4)+sk3
             u(n-ninc,5)=u(n-ninc,5)+sk4
          enddo
       enddo
!!$OMP END DO
!
!!$OMP DO
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k1)
          ind2 = indc(i2m1,j,k1)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             sk0= 2.*(v(n-ninc,2)*sn(m,kdir,1) &
                  +v(n-ninc,3)*sn(m,kdir,2) &
                  +v(n-ninc,4)*sn(m,kdir,3))
             sk1= 2.*(fxx(m-ninc)*sn(m,kdir,1) &
                  +fxy(m-ninc)*sn(m,kdir,2) &
                  +fxz(m-ninc)*sn(m,kdir,3))
             sk2= 2.*(fxy(m-ninc)*sn(m,kdir,1) &
                  +fyy(m-ninc)*sn(m,kdir,2) &
                  +fyz(m-ninc)*sn(m,kdir,3))
             sk3= 2.*(fxz(m-ninc)*sn(m,kdir,1) &
                  +fyz(m-ninc)*sn(m,kdir,2) &
                  +fzz(m-ninc)*sn(m,kdir,3))
             sk4= 2.*(fex(m-ninc)*sn(m,kdir,1) &
                  +fey(m-ninc)*sn(m,kdir,2) &
                  +fez(m-ninc)*sn(m,kdir,3))
             u(n,1)=u(n,1)-sk0
             u(n,2)=u(n,2)-sk1
             u(n,3)=u(n,3)-sk2
             u(n,4)=u(n,4)-sk3
             u(n,5)=u(n,5)-sk4
          enddo
       enddo
!!$OMP END DO
!
!!$OMP DO
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k2)
          ind2 = indc(i2m1,j,k2)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             sk0= 2.*(v(n,2)*sn(m,kdir,1) &
                  +v(n,3)*sn(m,kdir,2) &
                  +v(n,4)*sn(m,kdir,3))
             sk1= 2.*(fxx(m)*sn(m,kdir,1) &
                  +fxy(m)*sn(m,kdir,2) &
                  +fxz(m)*sn(m,kdir,3))
             sk2= 2.*(fxy(m)*sn(m,kdir,1) &
                  +fyy(m)*sn(m,kdir,2) &
                  +fyz(m)*sn(m,kdir,3))
             sk3= 2.*(fxz(m)*sn(m,kdir,1) &
                  +fyz(m)*sn(m,kdir,2) &
                  +fzz(m)*sn(m,kdir,3))
             sk4= 2.*(fex(m)*sn(m,kdir,1) &
                  +fey(m)*sn(m,kdir,2) &
                  +fez(m)*sn(m,kdir,3))
             u(n-ninc,1)=u(n-ninc,1)+sk0
             u(n-ninc,2)=u(n-ninc,2)+sk1
             u(n-ninc,3)=u(n-ninc,3)+sk2
             u(n-ninc,4)=u(n-ninc,4)+sk3
             u(n-ninc,5)=u(n-ninc,5)+sk4
          enddo
       enddo
!!$OMP END DO
!
    endif
!
!------normalisation et ajout de la dissipation artificielle------
!
!!$OMP DO collapse(2)
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1  ,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             u(n,1)=0.5*u(n,1)-d(n,1)
             u(n,2)=0.5*u(n,2)-d(n,2)
             u(n,3)=0.5*u(n,3)-d(n,3)
             u(n,4)=0.5*u(n,4)-d(n,4)
             u(n,5)=0.5*u(n,5)-d(n,5)
          enddo
       enddo
    enddo
!!$OMP END DO
!
!      write(6,'("===>sch_jam4: ecriture increment expli")')
!      k=1
!      i=40
!      do j=j1,j2m1
!       n=indc(i,j,k)
!       m=n-n0c
!       write(6,'(i4,i6,4(1pe12.4))')
!     &    j,n,u(n,1),u(n,2),u(n,3),u(n,5)
!      enddo
!
!-----calcul de la 'forcing function'---------------------------
!
    if(ityprk.ne.0) then
!!$OMP DO collapse(2)
       do k=k1,k2m1
          do j=j1,j2m1
             ind1=indc(i1  ,j,k)
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
!!$OMP END DO
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
  end subroutine sch_jameson3
end module mod_sch_jameson3
