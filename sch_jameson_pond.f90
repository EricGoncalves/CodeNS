module mod_sch_jameson_pond
implicit none
contains
      subroutine sch_jameson_pond( &
                 lm,ityprk, &
                 u,v,d,ff, &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 equat, &
                 sn,lgsnlt, &
                 fxx,fyy,fzz,fxy,fxz,fyz,fex,fey,fez, &
                 ps, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DA  DATE_C : mars 2002 - Eric Goncalves / Sinumef
!
!     ACT
!_A    Calcul des bilans de flux physiques pour chaque cellule.
!_A    Schema de Jameson avec ponderation pour maillage irregulier.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use proprieteflu
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
double precision :: d
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
double precision :: fxx
double precision :: fyy
double precision :: fzz
double precision :: fxy
double precision :: fxz
double precision :: fyz
double precision :: fex
double precision :: fey
double precision :: fez
double precision :: ps
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
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
integer :: nci
integer :: ncj
integer :: nck
integer :: nid
integer :: nijd
integer :: ninc
integer :: njd
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
!
!-----------------------------------------------------------------------
!
      character(len=7 ) :: equat
      integer isortie
      dimension v(ip11,ip60),d(ip11,ip60),u(ip11,ip60),ff(ip11,ip60)
      dimension toxx(ip12),toxy(ip12),toxz(ip12), &
                toyy(ip12),toyz(ip12),tozz(ip12), &
                qcx(ip12),qcy(ip12),qcz(ip12)
      dimension sn(lgsnlt,nind,ndir)
      dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21),ps(ip11)
      dimension fxx(ip00),fyy(ip00),fzz(ip00),fxy(ip00),fxz(ip00), &
                fyz(ip00),fex(ip00),fey(ip00),fez(ip00)
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
!-----initalisation--------------------------------
!
      ind1 = indc(i1m1,j1m1,k1m1)
      ind2 = indc(i2  ,j2  ,k2  )
      do n=ind1,ind2
       u(n,1)=0.
       u(n,2)=0.
       u(n,3)=0.
       u(n,4)=0.
       u(n,5)=0.
      enddo
!
!-----calcul des densites de flux convectifs et visqueux--------------
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
!*********************************************************************
!      calcul des flux numeriques par direction
!*********************************************************************
!
!------direction i----------------------------------------------
!
      kdir=1
      ninc=nci
!
      do k=k1,k2m1
       do j=j1,j2m1
        ind1 = indc(i1p1,j,k)
        ind2 = indc(i2m1,j,k)
!DEC$ IVDEP
        do n=ind1,ind2
         m=n-n0c
         si0=(cmui1(m)*v(n,2)+cmui2(m)*v(n-ninc,2))*sn(m,kdir,1) &
            +(cmui1(m)*v(n,3)+cmui2(m)*v(n-ninc,3))*sn(m,kdir,2) &
            +(cmui1(m)*v(n,4)+cmui2(m)*v(n-ninc,4))*sn(m,kdir,3)
         si1=(cmui1(m)*fxx(m)+cmui2(m)*fxx(m-ninc))*sn(m,kdir,1) &
            +(cmui1(m)*fxy(m)+cmui2(m)*fxy(m-ninc))*sn(m,kdir,2) &
            +(cmui1(m)*fxz(m)+cmui2(m)*fxz(m-ninc))*sn(m,kdir,3)
         si2=(cmui1(m)*fxy(m)+cmui2(m)*fxy(m-ninc))*sn(m,kdir,1) &
            +(cmui1(m)*fyy(m)+cmui2(m)*fyy(m-ninc))*sn(m,kdir,2) &
            +(cmui1(m)*fyz(m)+cmui2(m)*fyz(m-ninc))*sn(m,kdir,3)
         si3=(cmui1(m)*fxz(m)+cmui2(m)*fxz(m-ninc))*sn(m,kdir,1) &
            +(cmui1(m)*fyz(m)+cmui2(m)*fyz(m-ninc))*sn(m,kdir,2) &
            +(cmui1(m)*fzz(m)+cmui2(m)*fzz(m-ninc))*sn(m,kdir,3)
         si4=(cmui1(m)*fex(m)+cmui2(m)*fex(m-ninc))*sn(m,kdir,1) &
            +(cmui1(m)*fey(m)+cmui2(m)*fey(m-ninc))*sn(m,kdir,2) &
            +(cmui1(m)*fez(m)+cmui2(m)*fez(m-ninc))*sn(m,kdir,3)
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
       ind1 = indc(i1,j1  ,k)
       ind2 = indc(i1,j2m1,k)
!DEC$ IVDEP
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
!DEC$ IVDEP
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
!------direction j----------------------------------------------
!
      kdir=2
      ninc=ncj
!
      do k=k1,k2m1
       do j=j1p1,j2m1
        ind1 = indc(i1  ,j,k)
        ind2 = indc(i2m1,j,k)
!DEC$ IVDEP
        do n=ind1,ind2
         m=n-n0c
         sj0=(cmuj1(m)*v(n,2)+cmuj2(m)*v(n-ninc,2))*sn(m,kdir,1) &
            +(cmuj1(m)*v(n,3)+cmuj2(m)*v(n-ninc,3))*sn(m,kdir,2) &
            +(cmuj1(m)*v(n,4)+cmuj2(m)*v(n-ninc,4))*sn(m,kdir,3)
         sj1=(cmuj1(m)*fxx(m)+cmuj2(m)*fxx(m-ninc))*sn(m,kdir,1) &
            +(cmuj1(m)*fxy(m)+cmuj2(m)*fxy(m-ninc))*sn(m,kdir,2) &
            +(cmuj1(m)*fxz(m)+cmuj2(m)*fxz(m-ninc))*sn(m,kdir,3)
         sj2=(cmuj1(m)*fxy(m)+cmuj2(m)*fxy(m-ninc))*sn(m,kdir,1) &
            +(cmuj1(m)*fyy(m)+cmuj2(m)*fyy(m-ninc))*sn(m,kdir,2) &
            +(cmuj1(m)*fyz(m)+cmuj2(m)*fyz(m-ninc))*sn(m,kdir,3)
         sj3=(cmuj1(m)*fxz(m)+cmuj2(m)*fxz(m-ninc))*sn(m,kdir,1) &
            +(cmuj1(m)*fyz(m)+cmuj2(m)*fyz(m-ninc))*sn(m,kdir,2) &
            +(cmuj1(m)*fzz(m)+cmuj2(m)*fzz(m-ninc))*sn(m,kdir,3)
         sj4=(cmuj1(m)*fex(m)+cmuj2(m)*fex(m-ninc))*sn(m,kdir,1) &
            +(cmuj1(m)*fey(m)+cmuj2(m)*fey(m-ninc))*sn(m,kdir,2) &
            +(cmuj1(m)*fez(m)+cmuj2(m)*fez(m-ninc))*sn(m,kdir,3)
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
       ind1 = indc(i1  ,j1,k)
       ind2 = indc(i2m1,j1,k)
!DEC$ IVDEP
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
!DEC$ IVDEP
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
      do k=k1p1,k2m1
       do j=j1,j2m1
        ind1 = indc(i1  ,j,k)
        ind2 = indc(i2m1,j,k)
!DEC$ IVDEP
        do n=ind1,ind2
         m=n-n0c
         sk0=(cmuk1(m)*v(n,2)+cmuk2(m)*v(n-ninc,2))*sn(m,kdir,1) &
            +(cmuk1(m)*v(n,3)+cmuk2(m)*v(n-ninc,3))*sn(m,kdir,2) &
            +(cmuk1(m)*v(n,4)+cmuk2(m)*v(n-ninc,4))*sn(m,kdir,3)
         sk1=(cmuk1(m)*fxx(m)+cmuk2(m)*fxx(m-ninc))*sn(m,kdir,1) &
            +(cmuk1(m)*fxy(m)+cmuk2(m)*fxy(m-ninc))*sn(m,kdir,2) &
            +(cmuk1(m)*fxz(m)+cmuk2(m)*fxz(m-ninc))*sn(m,kdir,3)
         sk2=(cmuk1(m)*fxy(m)+cmuk2(m)*fxy(m-ninc))*sn(m,kdir,1) &
            +(cmuk1(m)*fyy(m)+cmuk2(m)*fyy(m-ninc))*sn(m,kdir,2) &
            +(cmuk1(m)*fyz(m)+cmuk2(m)*fyz(m-ninc))*sn(m,kdir,3)
         sk3=(cmuk1(m)*fxz(m)+cmuk2(m)*fxz(m-ninc))*sn(m,kdir,1) &
            +(cmuk1(m)*fyz(m)+cmuk2(m)*fyz(m-ninc))*sn(m,kdir,2) &
            +(cmuk1(m)*fzz(m)+cmuk2(m)*fzz(m-ninc))*sn(m,kdir,3)
         sk4=(cmuk1(m)*fex(m)+cmuk2(m)*fex(m-ninc))*sn(m,kdir,1) &
            +(cmuk1(m)*fey(m)+cmuk2(m)*fey(m-ninc))*sn(m,kdir,2) &
            +(cmuk1(m)*fez(m)+cmuk2(m)*fez(m-ninc))*sn(m,kdir,3)
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
       ind1 = indc(i1  ,j,k1)
       ind2 = indc(i2m1,j,k1)
!DEC$ IVDEP
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
!DEC$ IVDEP
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
         ind1=indc(i1,j,k)
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
      isortie=0
      if(isortie.eq.1) then
      write(6,'("===>sch_jam3_pond: ecriture increment expli")')
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
