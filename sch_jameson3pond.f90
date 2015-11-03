module mod_sch_jameson3pond
  implicit none
contains
  subroutine sch_jameson3pond( &
       lm,ityprk, &
       u,v,d,ff, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       equat, &
       sn,lgsnlt, &
       fxx,fyy,fzz,fxy,fxz,fyz,fex,fey,fez, &
       ps, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
       cvi,cvj,cvk)
!
!*****************************************************************
!
!_DA  DATE_C : avril 2002 - Eric Goncalves / Sinumef
!
!     ACT
!_A    Calcul des bilans de flux physiques pour chaque cellule.
!_A    Schema de Jameson avec correction de l'erreur dispersive.
!_A    Ordre 3 en cartesien.
!
!*******************************************************************
!-----parameters figes----------------------------------------------
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
    double precision ::                   a1,                  a2,                  a3,                  a4,         cmui1(ip21)
    double precision ::          cmui2(ip21),         cmuj1(ip21),         cmuj2(ip21),         cmuk1(ip21),         cmuk2(ip21)
    double precision ::            cvi(ip21),           cvj(ip21),           cvk(ip21),        d(ip11,ip60),           fex(ip00)
    double precision ::            fey(ip00),           fez(ip00),       ff(ip11,ip60),           fxx(ip00),           fxy(ip00)
    double precision ::            fxz(ip00),           fyy(ip00),           fyz(ip00),           fzz(ip00),            ps(ip11)
    double precision ::            qcx(ip12),           qcy(ip12),           qcz(ip12),                 si0,                 si1
    double precision ::                  si2,                 si3,                 si4,                 sj0,                 sj1
    double precision ::                  sj2,                 sj3,                 sj4,                 sk0,                 sk1
    double precision ::                  sk2,                 sk3,                 sk4,sn(lgsnlt,nind,ndir),          toxx(ip12)
    double precision ::           toxy(ip12),          toxz(ip12),          toyy(ip12),          toyz(ip12),          tozz(ip12)
    double precision ::         u(ip11,ip60),        v(ip11,ip60)
!
!-----------------------------------------------------------------
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
!-----calcul des densites de flux visqueuses---------------------
!
    if(equat(3:5).eq.'2dk') then
       ind1 = indc(i1m1,j1m1,k1  )
       ind2 = indc(i2  ,j2  ,k2m1)
    elseif(equat(3:4).eq.'3d') then
       ind1 = indc(i1m1,j1m1,k1m1)
       ind2 = indc(i2  ,j2  ,k2  )
    endif
!
    if (equat(1:2).eq.'ns') then
       do n=ind1,ind2
          m=n-n0c
          u(n,1)=0.
          u(n,2)=0.
          u(n,3)=0.
          u(n,4)=0.
          u(n,5)=0.
          fxx(m)=v(n,2)*(v(n,2)/v(n,1))+ps(n)-pinfl-toxx(n)
          fxy(m)=v(n,3)*(v(n,2)/v(n,1))  -toxy(n)
          fxz(m)=v(n,4)*(v(n,2)/v(n,1))  -toxz(n)
          fyy(m)=v(n,3)*(v(n,3)/v(n,1))+ps(n)-pinfl-toyy(n)
          fyz(m)=v(n,4)*(v(n,3)/v(n,1))  -toyz(n)
          fzz(m)=v(n,4)*(v(n,4)/v(n,1))+ps(n)-pinfl-tozz(n)
          fex(m)=((v(n,5)+ps(n)-pinfl-toxx(n))*v(n,2) &
               -toxy(n)*v(n,3)-toxz(n)*v(n,4))/v(n,1)-qcx(n)
          fey(m)=((v(n,5)+ps(n)-pinfl-toyy(n))*v(n,3) &
               -toxy(n)*v(n,2)-toyz(n)*v(n,4))/v(n,1)-qcy(n)
          fez(m)=((v(n,5)+ps(n)-pinfl-tozz(n))*v(n,4) &
               -toxz(n)*v(n,2)-toyz(n)*v(n,3))/v(n,1)-qcz(n)
       enddo
    else
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
    endif
!
!**************************************************************
!      calcul des flux numeriques par direction
!**********************************************************
!
!------direction i-------------------------------------
!
    kdir=1
    ninc=nci
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1+2,j,k)
          ind2 = indc(i2-2,j,k)
          do n=ind1,ind2
             m=n-n0c
             n2=n-2*ninc
             m2=m-2*ninc
             a1=-0.25*cmui1(n)*cmui2(n)*cvi(n)**2/(cvi(n+ninc)*(cvi(n)+cvi(n+ninc)))
             a2=cmui1(n)+0.25*cmui1(n)*cmui2(n)*cvi(n)*(1./cvi(n+ninc) &
                  -1./(cvi(n-ninc)+cvi(n)))
             a3=cmui2(n)+0.25*cmui1(n)*cmui2(n)*cvi(n)*(1./cvi(n-ninc) &
                  -1./(cvi(n)+cvi(n+ninc)))
             a4=-0.25*cmui1(n)*cmui2(n)*cvi(n)**2/(cvi(n-ninc)*(cvi(n-ninc)+cvi(n)))
!
             si0=(a1*v(n+ninc,2)+a2*v(n,2)+a3*v(n-ninc,2)+a4*v(n2,2))*sn(m,kdir,1) &
                +(a1*v(n+ninc,3)+a2*v(n,3)+a3*v(n-ninc,3)+a4*v(n2,3))*sn(m,kdir,2) &
                +(a1*v(n+ninc,4)+a2*v(n,4)+a3*v(n-ninc,4)+a4*v(n2,4))*sn(m,kdir,3)
             si1=(a1*fxx(m+ninc)+a2*fxx(m)+a3*fxx(m-ninc)+a4*fxx(m2))*sn(m,kdir,1) &
                +(a1*fxy(m+ninc)+a2*fxy(m)+a3*fxy(m-ninc)+a4*fxy(m2))*sn(m,kdir,2) &
                +(a1*fxz(m+ninc)+a2*fxz(m)+a3*fxz(m-ninc)+a4*fxz(m2))*sn(m,kdir,3)
             si2=(a1*fxy(m+ninc)+a2*fxy(m)+a3*fxy(m-ninc)+a4*fxy(m2))*sn(m,kdir,1) &
                +(a1*fyy(m+ninc)+a2*fyy(m)+a3*fyy(m-ninc)+a4*fyy(m2))*sn(m,kdir,2) &
                +(a1*fyz(m+ninc)+a2*fyz(m)+a3*fyz(m-ninc)+a4*fyz(m2))*sn(m,kdir,3)
             si3=(a1*fxz(m+ninc)+a2*fxz(m)+a3*fxz(m-ninc)+a4*fxz(m2))*sn(m,kdir,1) &
                +(a1*fyz(m+ninc)+a2*fyz(m)+a3*fyz(m-ninc)+a4*fyz(m2))*sn(m,kdir,2) &
                +(a1*fzz(m+ninc)+a2*fzz(m)+a3*fzz(m-ninc)+a4*fzz(m2))*sn(m,kdir,3)
             si4=(a1*fex(m+ninc)+a2*fex(m)+a3*fex(m-ninc)+a4*fex(m2))*sn(m,kdir,1) &
                +(a1*fey(m+ninc)+a2*fey(m)+a3*fey(m-ninc)+a4*fey(m2))*sn(m,kdir,2) &
                +(a1*fez(m+ninc)+a2*fez(m)+a3*fez(m-ninc)+a4*fez(m2))*sn(m,kdir,3)
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
!
    do k=k1,k2m1
       ind1 = indc(i1p1,j1  ,k)
       ind2 = indc(i1p1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          si0=(cmui1(n)*v(n,2)+cmui2(n)*v(n-ninc,2))*sn(m,kdir,1) &
             +(cmui1(n)*v(n,3)+cmui2(n)*v(n-ninc,3))*sn(m,kdir,2) &
             +(cmui1(n)*v(n,4)+cmui2(n)*v(n-ninc,4))*sn(m,kdir,3)
          si1=(cmui1(n)*fxx(m)+cmui2(n)*fxx(m-ninc))*sn(m,kdir,1) &
             +(cmui1(n)*fxy(m)+cmui2(n)*fxy(m-ninc))*sn(m,kdir,2) &
             +(cmui1(n)*fxz(m)+cmui2(n)*fxz(m-ninc))*sn(m,kdir,3)
          si2=(cmui1(n)*fxy(m)+cmui2(n)*fxy(m-ninc))*sn(m,kdir,1) &
             +(cmui1(n)*fyy(m)+cmui2(n)*fyy(m-ninc))*sn(m,kdir,2) &
             +(cmui1(n)*fyz(m)+cmui2(n)*fyz(m-ninc))*sn(m,kdir,3)
          si3=(cmui1(n)*fxz(m)+cmui2(n)*fxz(m-ninc))*sn(m,kdir,1) &
             +(cmui1(n)*fyz(m)+cmui2(n)*fyz(m-ninc))*sn(m,kdir,2) &
             +(cmui1(n)*fzz(m)+cmui2(n)*fzz(m-ninc))*sn(m,kdir,3)
          si4=(cmui1(n)*fex(m)+cmui2(n)*fex(m-ninc))*sn(m,kdir,1) &
             +(cmui1(n)*fey(m)+cmui2(n)*fey(m-ninc))*sn(m,kdir,2) &
             +(cmui1(n)*fez(m)+cmui2(n)*fez(m-ninc))*sn(m,kdir,3)
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
!
    do k=k1,k2m1
       ind1 = indc(i2m1,j1  ,k)
       ind2 = indc(i2m1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          si0=(cmui1(n)*v(n,2)+cmui2(n)*v(n-ninc,2))*sn(m,kdir,1) &
             +(cmui1(n)*v(n,3)+cmui2(n)*v(n-ninc,3))*sn(m,kdir,2) &
             +(cmui1(n)*v(n,4)+cmui2(n)*v(n-ninc,4))*sn(m,kdir,3)
          si1=(cmui1(n)*fxx(m)+cmui2(n)*fxx(m-ninc))*sn(m,kdir,1) &
             +(cmui1(n)*fxy(m)+cmui2(n)*fxy(m-ninc))*sn(m,kdir,2) &
             +(cmui1(n)*fxz(m)+cmui2(n)*fxz(m-ninc))*sn(m,kdir,3)
          si2=(cmui1(n)*fxy(m)+cmui2(n)*fxy(m-ninc))*sn(m,kdir,1) &
             +(cmui1(n)*fyy(m)+cmui2(n)*fyy(m-ninc))*sn(m,kdir,2) &
             +(cmui1(n)*fyz(m)+cmui2(n)*fyz(m-ninc))*sn(m,kdir,3)
          si3=(cmui1(n)*fxz(m)+cmui2(n)*fxz(m-ninc))*sn(m,kdir,1) &
             +(cmui1(n)*fyz(m)+cmui2(n)*fyz(m-ninc))*sn(m,kdir,2) &
             +(cmui1(n)*fzz(m)+cmui2(n)*fzz(m-ninc))*sn(m,kdir,3)
          si4=(cmui1(n)*fex(m)+cmui2(n)*fex(m-ninc))*sn(m,kdir,1) &
             +(cmui1(n)*fey(m)+cmui2(n)*fey(m-ninc))*sn(m,kdir,2) &
             +(cmui1(n)*fez(m)+cmui2(n)*fez(m-ninc))*sn(m,kdir,3)
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
!
    do k=k1,k2m1
       ind1 = indc(i1,j1  ,k)
       ind2 = indc(i1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          si0= 2*v(n-ninc,2)*sn(m,kdir,1) &
              +2*v(n-ninc,3)*sn(m,kdir,2) &
              +2*v(n-ninc,4)*sn(m,kdir,3)
          si1= 2*fxx(m-ninc)*sn(m,kdir,1) &
              +2*fxy(m-ninc)*sn(m,kdir,2) &
              +2*fxz(m-ninc)*sn(m,kdir,3)
          si2= 2*fxy(m-ninc)*sn(m,kdir,1) &
              +2*fyy(m-ninc)*sn(m,kdir,2) &
              +2*fyz(m-ninc)*sn(m,kdir,3)
          si3= 2*fxz(m-ninc)*sn(m,kdir,1) &
              +2*fyz(m-ninc)*sn(m,kdir,2) &
              +2*fzz(m-ninc)*sn(m,kdir,3)
          si4= 2*fex(m-ninc)*sn(m,kdir,1) &
              +2*fey(m-ninc)*sn(m,kdir,2) &
              +2*fez(m-ninc)*sn(m,kdir,3)
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
          si0= 2*v(n,2)*sn(m,kdir,1) &
              +2*v(n,3)*sn(m,kdir,2) &
              +2*v(n,4)*sn(m,kdir,3)
          si1= 2*fxx(m)*sn(m,kdir,1) &
              +2*fxy(m)*sn(m,kdir,2) &
              +2*fxz(m)*sn(m,kdir,3)
          si2= 2*fxy(m)*sn(m,kdir,1) &
              +2*fyy(m)*sn(m,kdir,2) &
              +2*fyz(m)*sn(m,kdir,3)
          si3= 2*fxz(m)*sn(m,kdir,1) &
              +2*fyz(m)*sn(m,kdir,2) &
              +2*fzz(m)*sn(m,kdir,3)
          si4= 2*fex(m)*sn(m,kdir,1) &
              +2*fey(m)*sn(m,kdir,2) &
              +2*fez(m)*sn(m,kdir,3)
          u(n-ninc,1)=u(n-ninc,1)+si0
          u(n-ninc,2)=u(n-ninc,2)+si1
          u(n-ninc,3)=u(n-ninc,3)+si2
          u(n-ninc,4)=u(n-ninc,4)+si3
          u(n-ninc,5)=u(n-ninc,5)+si4
       enddo
    enddo
!
!------direction j----------------------------------------
!
    kdir=2
    ninc=ncj
!
    do k=k1,k2m1
       do j=j1+2,j2-2
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             n2=n-2*ninc
             m2=m-2*ninc
             a1=-0.25*cmuj1(n)*cmuj2(n)*cvj(n)**2/(cvj(n+ninc)*(cvj(n)+cvj(n+ninc)))
             a2=cmuj1(n)+0.25*cmuj1(n)*cmuj2(n)*cvj(n)*(1./cvj(n+ninc) &
                  -1./(cvj(n-ninc)+cvj(n)))
             a3=cmuj2(n)+0.25*cmuj1(n)*cmuj2(n)*cvj(n)*(1./cvj(n-ninc) &
                  -1./(cvj(n)+cvj(n+ninc)))
             a4=-0.25*cmuj1(n)*cmuj2(n)*cvj(n)**2/(cvj(n-ninc)*(cvj(n-ninc)+cvj(n)))
!
             sj0=(a1*v(n+ninc,2)+a2*v(n,2)+a3*v(n-ninc,2)+a4*v(n2,2))*sn(m,kdir,1) &
                +(a1*v(n+ninc,3)+a2*v(n,3)+a3*v(n-ninc,3)+a4*v(n2,3))*sn(m,kdir,2) &
                +(a1*v(n+ninc,4)+a2*v(n,4)+a3*v(n-ninc,4)+a4*v(n2,4))*sn(m,kdir,3)
             sj1=(a1*fxx(m+ninc)+a2*fxx(m)+a3*fxx(m-ninc)+a4*fxx(m2))*sn(m,kdir,1) &
                +(a1*fxy(m+ninc)+a2*fxy(m)+a3*fxy(m-ninc)+a4*fxy(m2))*sn(m,kdir,2) &
                +(a1*fxz(m+ninc)+a2*fxz(m)+a3*fxz(m-ninc)+a4*fxz(m2))*sn(m,kdir,3)
             sj2=(a1*fxy(m+ninc)+a2*fxy(m)+a3*fxy(m-ninc)+a4*fxy(m2))*sn(m,kdir,1) &
                +(a1*fyy(m+ninc)+a2*fyy(m)+a3*fyy(m-ninc)+a4*fyy(m2))*sn(m,kdir,2) &
                +(a1*fyz(m+ninc)+a2*fyz(m)+a3*fyz(m-ninc)+a4*fyz(m2))*sn(m,kdir,3)
             sj3=(a1*fxz(m+ninc)+a2*fxz(m)+a3*fxz(m-ninc)+a4*fxz(m2))*sn(m,kdir,1) &
                +(a1*fyz(m+ninc)+a2*fyz(m)+a3*fyz(m-ninc)+a4*fyz(m2))*sn(m,kdir,2) &
                +(a1*fzz(m+ninc)+a2*fzz(m)+a3*fzz(m-ninc)+a4*fzz(m2))*sn(m,kdir,3)
             sj4=(a1*fex(m+ninc)+a2*fex(m)+a3*fex(m-ninc)+a4*fex(m2))*sn(m,kdir,1) &
                +(a1*fey(m+ninc)+a2*fey(m)+a3*fey(m-ninc)+a4*fey(m2))*sn(m,kdir,2) &
                +(a1*fez(m+ninc)+a2*fez(m)+a3*fez(m-ninc)+a4*fez(m2))*sn(m,kdir,3)
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
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j1p1,k)
       ind2 = indc(i2m1,j1p1,k)
       do n=ind1,ind2
          m=n-n0c
          sj0=(cmuj1(n)*v(n,2)+cmuj2(n)*v(n-ninc,2))*sn(m,kdir,1) &
             +(cmuj1(n)*v(n,3)+cmuj2(n)*v(n-ninc,3))*sn(m,kdir,2) &
             +(cmuj1(n)*v(n,4)+cmuj2(n)*v(n-ninc,4))*sn(m,kdir,3)
          sj1=(cmuj1(n)*fxx(m)+cmuj2(n)*fxx(m-ninc))*sn(m,kdir,1) &
             +(cmuj1(n)*fxy(m)+cmuj2(n)*fxy(m-ninc))*sn(m,kdir,2) &
             +(cmuj1(n)*fxz(m)+cmuj2(n)*fxz(m-ninc))*sn(m,kdir,3)
          sj2=(cmuj1(n)*fxy(m)+cmuj2(n)*fxy(m-ninc))*sn(m,kdir,1) &
             +(cmuj1(n)*fyy(m)+cmuj2(n)*fyy(m-ninc))*sn(m,kdir,2) &
             +(cmuj1(n)*fyz(m)+cmuj2(n)*fyz(m-ninc))*sn(m,kdir,3)
          sj3=(cmuj1(n)*fxz(m)+cmuj2(n)*fxz(m-ninc))*sn(m,kdir,1) &
             +(cmuj1(n)*fyz(m)+cmuj2(n)*fyz(m-ninc))*sn(m,kdir,2) &
             +(cmuj1(n)*fzz(m)+cmuj2(n)*fzz(m-ninc))*sn(m,kdir,3)
          sj4=(cmuj1(n)*fex(m)+cmuj2(n)*fex(m-ninc))*sn(m,kdir,1) &
             +(cmuj1(n)*fey(m)+cmuj2(n)*fey(m-ninc))*sn(m,kdir,2) &
             +(cmuj1(n)*fez(m)+cmuj2(n)*fez(m-ninc))*sn(m,kdir,3)
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
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j2m1,k)
       ind2 = indc(i2m1,j2m1,k)
       do n=ind1,ind2
          m=n-n0c
          sj0=(cmuj1(n)*v(n,2)+cmuj2(n)*v(n-ninc,2))*sn(m,kdir,1) &
             +(cmuj1(n)*v(n,3)+cmuj2(n)*v(n-ninc,3))*sn(m,kdir,2) &
             +(cmuj1(n)*v(n,4)+cmuj2(n)*v(n-ninc,4))*sn(m,kdir,3)
          sj1=(cmuj1(n)*fxx(m)+cmuj2(n)*fxx(m-ninc))*sn(m,kdir,1) &
             +(cmuj1(n)*fxy(m)+cmuj2(n)*fxy(m-ninc))*sn(m,kdir,2) &
             +(cmuj1(n)*fxz(m)+cmuj2(n)*fxz(m-ninc))*sn(m,kdir,3)
          sj2=(cmuj1(n)*fxy(m)+cmuj2(n)*fxy(m-ninc))*sn(m,kdir,1) &
             +(cmuj1(n)*fyy(m)+cmuj2(n)*fyy(m-ninc))*sn(m,kdir,2) &
             +(cmuj1(n)*fyz(m)+cmuj2(n)*fyz(m-ninc))*sn(m,kdir,3)
          sj3=(cmuj1(n)*fxz(m)+cmuj2(n)*fxz(m-ninc))*sn(m,kdir,1) &
             +(cmuj1(n)*fyz(m)+cmuj2(n)*fyz(m-ninc))*sn(m,kdir,2) &
             +(cmuj1(n)*fzz(m)+cmuj2(n)*fzz(m-ninc))*sn(m,kdir,3)
          sj4=(cmuj1(n)*fex(m)+cmuj2(n)*fex(m-ninc))*sn(m,kdir,1) &
             +(cmuj1(n)*fey(m)+cmuj2(n)*fey(m-ninc))*sn(m,kdir,2) &
             +(cmuj1(n)*fez(m)+cmuj2(n)*fez(m-ninc))*sn(m,kdir,3)
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
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j1,k)
       ind2 = indc(i2m1,j1,k)
       do n=ind1,ind2
          m=n-n0c
          sj0= 2*v(n-ninc,2)*sn(m,kdir,1) &
              +2*v(n-ninc,3)*sn(m,kdir,2) &
              +2*v(n-ninc,4)*sn(m,kdir,3)
          sj1= 2*fxx(m-ninc)*sn(m,kdir,1) &
              +2*fxy(m-ninc)*sn(m,kdir,2) &
              +2*fxz(m-ninc)*sn(m,kdir,3)
          sj2= 2*fxy(m-ninc)*sn(m,kdir,1) &
              +2*fyy(m-ninc)*sn(m,kdir,2) &
              +2*fyz(m-ninc)*sn(m,kdir,3)
          sj3= 2*fxz(m-ninc)*sn(m,kdir,1) &
              +2*fyz(m-ninc)*sn(m,kdir,2) &
              +2*fzz(m-ninc)*sn(m,kdir,3)
          sj4= 2*fex(m-ninc)*sn(m,kdir,1) &
              +2*fey(m-ninc)*sn(m,kdir,2) &
              +2*fez(m-ninc)*sn(m,kdir,3)
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
          sj0= 2*v(n,2)*sn(m,kdir,1) &
              +2*v(n,3)*sn(m,kdir,2) &
              +2*v(n,4)*sn(m,kdir,3)
          sj1= 2*fxx(m)*sn(m,kdir,1) &
              +2*fxy(m)*sn(m,kdir,2) &
              +2*fxz(m)*sn(m,kdir,3)
          sj2= 2*fxy(m)*sn(m,kdir,1) &
              +2*fyy(m)*sn(m,kdir,2) &
              +2*fyz(m)*sn(m,kdir,3)
          sj3= 2*fxz(m)*sn(m,kdir,1) &
              +2*fyz(m)*sn(m,kdir,2) &
              +2*fzz(m)*sn(m,kdir,3)
          sj4= 2*fex(m)*sn(m,kdir,1) &
              +2*fey(m)*sn(m,kdir,2) &
              +2*fez(m)*sn(m,kdir,3)
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
       do k=k1+2,k2-2
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                n2=n-2*ninc
                m2=m-2*ninc
                a1=-0.25*cmuk1(n)*cmuk2(n)*cvk(n)**2/(cvk(n+ninc)*(cvk(n)+cvk(n+ninc)))
                a2=cmuk1(n)+0.25*cmuk1(n)*cmuk2(n)*cvk(n)*(1./cvk(n+ninc) &
                     -1./(cvk(n-ninc)+cvk(n)))
                a3=cmuk2(n)+0.25*cmuk1(n)*cmuk2(n)*cvk(n)*(1./cvk(n-ninc) &
                     -1./(cvk(n)+cvk(n+ninc)))
                a4=-0.25*cmuk1(n)*cmuk2(n)*cvk(n)**2/(cvk(n-ninc)*(cvk(n-ninc)+cvk(n)))
!
                sk0=(a1*v(n+ninc,2)+a2*v(n,2)+a3*v(n-ninc,2)+a4*v(n2,2))*sn(m,kdir,1) &
                   +(a1*v(n+ninc,3)+a2*v(n,3)+a3*v(n-ninc,3)+a4*v(n2,3))*sn(m,kdir,2) &
                   +(a1*v(n+ninc,4)+a2*v(n,4)+a3*v(n-ninc,4)+a4*v(n2,4))*sn(m,kdir,3)
                sk1=(a1*fxx(m+ninc)+a2*fxx(m)+a3*fxx(m-ninc)+a4*fxx(m2))*sn(m,kdir,1) &
                   +(a1*fxy(m+ninc)+a2*fxy(m)+a3*fxy(m-ninc)+a4*fxy(m2))*sn(m,kdir,2) &
                   +(a1*fxz(m+ninc)+a2*fxz(m)+a3*fxz(m-ninc)+a4*fxz(m2))*sn(m,kdir,3)
                sk2=(a1*fxy(m+ninc)+a2*fxy(m)+a3*fxy(m-ninc)+a4*fxy(m2))*sn(m,kdir,1) &
                   +(a1*fyy(m+ninc)+a2*fyy(m)+a3*fyy(m-ninc)+a4*fyy(m2))*sn(m,kdir,2) &
                   +(a1*fyz(m+ninc)+a2*fyz(m)+a3*fyz(m-ninc)+a4*fyz(m2))*sn(m,kdir,3)
                sk3=(a1*fxz(m+ninc)+a2*fxz(m)+a3*fxz(m-ninc)+a4*fxz(m2))*sn(m,kdir,1) &
                   +(a1*fyz(m+ninc)+a2*fyz(m)+a3*fyz(m-ninc)+a4*fyz(m2))*sn(m,kdir,2) &
                   +(a1*fzz(m+ninc)+a2*fzz(m)+a3*fzz(m-ninc)+a4*fzz(m2))*sn(m,kdir,3)
                sk4=(a1*fex(m+ninc)+a2*fex(m)+a3*fex(m-ninc)+a4*fex(m2))*sn(m,kdir,1) &
                   +(a1*fey(m+ninc)+a2*fey(m)+a3*fey(m-ninc)+a4*fey(m2))*sn(m,kdir,2) &
                   +(a1*fez(m+ninc)+a2*fez(m)+a3*fez(m-ninc)+a4*fez(m2))*sn(m,kdir,3)
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
!
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k1p1)
          ind2 = indc(i2m1,j,k1p1)
          do n=ind1,ind2
             m=n-n0c
             sk0=(cmuk1(n)*v(n,2)+cmuk2(n)*v(n-ninc,2))*sn(m,kdir,1) &
                +(cmuk1(n)*v(n,3)+cmuk2(n)*v(n-ninc,3))*sn(m,kdir,2) &
                +(cmuk1(n)*v(n,4)+cmuk2(n)*v(n-ninc,4))*sn(m,kdir,3)
             sk1=(cmuk1(n)*fxx(m)+cmuk2(n)*fxx(m-ninc))*sn(m,kdir,1) &
                +(cmuk1(n)*fxy(m)+cmuk2(n)*fxy(m-ninc))*sn(m,kdir,2) &
                +(cmuk1(n)*fxz(m)+cmuk2(n)*fxz(m-ninc))*sn(m,kdir,3)
             sk2=(cmuk1(n)*fxy(m)+cmuk2(n)*fxy(m-ninc))*sn(m,kdir,1) &
                +(cmuk1(n)*fyy(m)+cmuk2(n)*fyy(m-ninc))*sn(m,kdir,2) &
                +(cmuk1(n)*fyz(m)+cmuk2(n)*fyz(m-ninc))*sn(m,kdir,3)
             sk3=(cmuk1(n)*fxz(m)+cmuk2(n)*fxz(m-ninc))*sn(m,kdir,1) &
                +(cmuk1(n)*fyz(m)+cmuk2(n)*fyz(m-ninc))*sn(m,kdir,2) &
                +(cmuk1(n)*fzz(m)+cmuk2(n)*fzz(m-ninc))*sn(m,kdir,3)
             sk4=(cmuk1(n)*fex(m)+cmuk2(n)*fex(m-ninc))*sn(m,kdir,1) &
                +(cmuk1(n)*fey(m)+cmuk2(n)*fey(m-ninc))*sn(m,kdir,2) &
                +(cmuk1(n)*fez(m)+cmuk2(n)*fez(m-ninc))*sn(m,kdir,3)
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
!
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k2m1)
          ind2 = indc(i2m1,j,k2m1)
          do n=ind1,ind2
             m=n-n0c
             sk0=(cmuk1(n)*v(n,2)+cmuk2(n)*v(n-ninc,2))*sn(m,kdir,1) &
                +(cmuk1(n)*v(n,3)+cmuk2(n)*v(n-ninc,3))*sn(m,kdir,2) &
                +(cmuk1(n)*v(n,4)+cmuk2(n)*v(n-ninc,4))*sn(m,kdir,3)
             sk1=(cmuk1(n)*fxx(m)+cmuk2(n)*fxx(m-ninc))*sn(m,kdir,1) &
                +(cmuk1(n)*fxy(m)+cmuk2(n)*fxy(m-ninc))*sn(m,kdir,2) &
                +(cmuk1(n)*fxz(m)+cmuk2(n)*fxz(m-ninc))*sn(m,kdir,3)
             sk2=(cmuk1(n)*fxy(m)+cmuk2(n)*fxy(m-ninc))*sn(m,kdir,1) &
                +(cmuk1(n)*fyy(m)+cmuk2(n)*fyy(m-ninc))*sn(m,kdir,2) &
                +(cmuk1(n)*fyz(m)+cmuk2(n)*fyz(m-ninc))*sn(m,kdir,3)
             sk3=(cmuk1(n)*fxz(m)+cmuk2(n)*fxz(m-ninc))*sn(m,kdir,1) &
                +(cmuk1(n)*fyz(m)+cmuk2(n)*fyz(m-ninc))*sn(m,kdir,2) &
                +(cmuk1(n)*fzz(m)+cmuk2(n)*fzz(m-ninc))*sn(m,kdir,3)
             sk4=(cmuk1(n)*fex(m)+cmuk2(n)*fex(m-ninc))*sn(m,kdir,1) &
                +(cmuk1(n)*fey(m)+cmuk2(n)*fey(m-ninc))*sn(m,kdir,2) &
                +(cmuk1(n)*fez(m)+cmuk2(n)*fez(m-ninc))*sn(m,kdir,3)
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
!
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k1)
          ind2 = indc(i2m1,j,k1)
          do n=ind1,ind2
             m=n-n0c
             sk0= 2*v(n-ninc,2)*sn(m,kdir,1) &
                  +2*v(n-ninc,3)*sn(m,kdir,2) &
                  +2*v(n-ninc,4)*sn(m,kdir,3)
             sk1= 2*fxx(m-ninc)*sn(m,kdir,1) &
                  +2*fxy(m-ninc)*sn(m,kdir,2) &
                  +2*fxz(m-ninc)*sn(m,kdir,3)
             sk2= 2*fxy(m-ninc)*sn(m,kdir,1) &
                  +2*fyy(m-ninc)*sn(m,kdir,2) &
                  +2*fyz(m-ninc)*sn(m,kdir,3)
             sk3= 2*fxz(m-ninc)*sn(m,kdir,1) &
                  +2*fyz(m-ninc)*sn(m,kdir,2) &
                  +2*fzz(m-ninc)*sn(m,kdir,3)
             sk4= 2*fex(m-ninc)*sn(m,kdir,1) &
                  +2*fey(m-ninc)*sn(m,kdir,2) &
                  +2*fez(m-ninc)*sn(m,kdir,3)
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
             sk0= 2*v(n,2)*sn(m,kdir,1) &
                  +2*v(n,3)*sn(m,kdir,2) &
                  +2*v(n,4)*sn(m,kdir,3)
             sk1= 2*fxx(m)*sn(m,kdir,1) &
                  +2*fxy(m)*sn(m,kdir,2) &
                  +2*fxz(m)*sn(m,kdir,3)
             sk2= 2*fxy(m)*sn(m,kdir,1) &
                  +2*fyy(m)*sn(m,kdir,2) &
                  +2*fyz(m)*sn(m,kdir,3)
             sk3= 2*fxz(m)*sn(m,kdir,1) &
                  +2*fyz(m)*sn(m,kdir,2) &
                  +2*fzz(m)*sn(m,kdir,3)
             sk4= 2*fex(m)*sn(m,kdir,1) &
                  +2*fey(m)*sn(m,kdir,2) &
                  +2*fez(m)*sn(m,kdir,3)
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
!------normalisation et ajout de la dissipation artificielle------
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1  ,j,k)
          ind2=indc(i2m1,j,k)
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
!
!      write(6,'("===>sch_jam4pond: ecriture increment expli")')
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
  end subroutine sch_jameson3pond
end module mod_sch_jameson3pond
