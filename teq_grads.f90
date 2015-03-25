module mod_teq_grads
  implicit none
contains
  subroutine teq_grads( &
       l, &
       equat, &
       sn,lgsnlt, &
       vol, &
       s, &
       dsdx,dsdy,dsdz, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DATE: decembre 2001 - Eric GONCALVES / SINUMEF
!
!     ACT     calcul du gradient de la vitesse
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use schemanum
    implicit none
  integer          ::      i,    i1,  i1m1,  i1p1,    i2
  integer          ::   i2m1,  i2p1,    id,  ind1,  ind2
  integer          ::      j,    j1,  j1m1,  j1p1,    j2
  integer          ::   j2m1,  j2p1,    jd,     k,    k1
  integer          ::   k1m1,  k1p1,    k2,  k2m1,  k2p1
  integer          ::     kd,  kdir,     l,lgsnlt,     m
  integer          ::      n,   n0c,   nci,   ncj,   nck
  integer          ::    nid,  nijd,  ninc,   njd
  double precision ::                   c0,         cmui1(ip21),         cmui2(ip21),         cmuj1(ip21),         cmuj2(ip21)
  double precision ::          cmuk1(ip21),         cmuk2(ip21),          dsdx(ip00),          dsdy(ip00),          dsdz(ip00)
  double precision ::                  eps,             s(ip00),                 si1,                 si2,                 si3
  double precision ::                  sj1,                 sj2,                 sj3,                 sk1,                 sk2
  double precision ::                  sk3,sn(lgsnlt,nind,ndir),                  ts,           vol(ip11),                vols
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
!
!



    eps=0.000001
!
    n0c=npc(l)
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd= nid*njd
!
    i1p1=i1+1
    j1p1=j1+1
    k1p1=k1+1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
    i2p1=i2+1
    j2p1=j2+1
    k2p1=k2+1
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
!
    nci=inc(1,0,0)
    ncj=inc(0,1,0)
    nck=inc(0,0,1)
!
!     initialisation
    ind1=indc(i1m1,j1m1,k1m1)-n0c
    ind2=indc(i2p1,j2p1,k2p1)-n0c
    do m=ind1,ind2
       dsdx(m)=0.
       dsdy(m)=0.
       dsdz(m)=0.
    enddo
!
!     test sur le schema numerique
    if(ischema.eq.1) then
!--------------------------------------------------------------------
!     direction k (a travers les facettes k=cste)
!--------------------------------------------------------------------
       if (equat(3:5).ne.'2dk') then
!
          kdir=3
          ninc=nck
!
          do k=k1p1,k2m1
             do j=j1,j2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   dsdx(m)=dsdx(m)-sn(m,kdir,1)*(s(m)+s(m-ninc))
                   dsdy(m)=dsdy(m)-sn(m,kdir,2)*(s(m)+s(m-ninc))
                   dsdz(m)=dsdz(m)-sn(m,kdir,3)*(s(m)+s(m-ninc))
                enddo
                do n=ind1,ind2
                   m=n-n0c
                   dsdx(m-ninc)=dsdx(m-ninc)+sn(m,kdir,1)*(s(m)+s(m-ninc))
                   dsdy(m-ninc)=dsdy(m-ninc)+sn(m,kdir,2)*(s(m)+s(m-ninc))
                   dsdz(m-ninc)=dsdz(m-ninc)+sn(m,kdir,3)*(s(m)+s(m-ninc))
                enddo
             enddo
          enddo
!
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k1)
             ind2 = indc(i2m1,j,k1)
             do n=ind1,ind2
                m=n-n0c
                dsdx(m)=dsdx(m)-sn(m,kdir,1)*2*s(m-ninc)
                dsdy(m)=dsdy(m)-sn(m,kdir,2)*2*s(m-ninc)
                dsdz(m)=dsdz(m)-sn(m,kdir,3)*2*s(m-ninc)
             enddo
          enddo
!
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k2)
             ind2 = indc(i2m1,j,k2)
             do n=ind1,ind2
                m=n-n0c
                dsdx(m-ninc)=dsdx(m-ninc)+sn(m,kdir,1)*2*s(m)
                dsdy(m-ninc)=dsdy(m-ninc)+sn(m,kdir,2)*2*s(m)
                dsdz(m-ninc)=dsdz(m-ninc)+sn(m,kdir,3)*2*s(m)
             enddo
          enddo
!
       endif
!--------------------------------------------------------------
!     direction j (a travers les facettes j=cste)
!--------------------------------------------------------------
       if (equat(3:5).ne.'2dj') then
!
          ninc=ncj
          kdir=2
!
          do j=j1p1,j2m1
             do k=k1,k2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   dsdx(m)=dsdx(m)-sn(m,kdir,1)*(s(m)+s(m-ninc))
                   dsdy(m)=dsdy(m)-sn(m,kdir,2)*(s(m)+s(m-ninc))
                   dsdz(m)=dsdz(m)-sn(m,kdir,3)*(s(m)+s(m-ninc))
                enddo
                do n=ind1,ind2
                   m=n-n0c
                   dsdx(m-ninc)=dsdx(m-ninc)+sn(m,kdir,1)*(s(m)+s(m-ninc))
                   dsdy(m-ninc)=dsdy(m-ninc)+sn(m,kdir,2)*(s(m)+s(m-ninc))
                   dsdz(m-ninc)=dsdz(m-ninc)+sn(m,kdir,3)*(s(m)+s(m-ninc))
                enddo
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i1  ,j1,k)
             ind2 = indc(i2m1,j1,k)
             do n=ind1,ind2
                m=n-n0c
                dsdx(m)=dsdx(m)-sn(m,kdir,1)*2*s(m-ninc)
                dsdy(m)=dsdy(m)-sn(m,kdir,2)*2*s(m-ninc)
                dsdz(m)=dsdz(m)-sn(m,kdir,3)*2*s(m-ninc)
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i1  ,j2,k)
             ind2 = indc(i2m1,j2,k)
             do n=ind1,ind2
                m=n-n0c
                dsdx(m-ninc)=dsdx(m-ninc)+sn(m,kdir,1)*2*s(m)
                dsdy(m-ninc)=dsdy(m-ninc)+sn(m,kdir,2)*2*s(m)
                dsdz(m-ninc)=dsdz(m-ninc)+sn(m,kdir,3)*2*s(m)
             enddo
          enddo
!
       endif
!---------------------------------------------------------------------
!     direction i (a travers les facettes i=cste)
!---------------------------------------------------------------------
       if (equat(3:5).ne.'2di') then
!
          ninc=nci
          kdir=1
!
          do k=k1,k2m1
             do j=j1,j2m1
                ind1 = indc(i1p1,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   dsdx(m)=dsdx(m)-sn(m,kdir,1)*(s(m)+s(m-ninc))
                   dsdy(m)=dsdy(m)-sn(m,kdir,2)*(s(m)+s(m-ninc))
                   dsdz(m)=dsdz(m)-sn(m,kdir,3)*(s(m)+s(m-ninc))
                enddo
                do n=ind1,ind2
                   m=n-n0c
                   dsdx(m-ninc)=dsdx(m-ninc)+sn(m,kdir,1)*(s(m)+s(m-ninc))
                   dsdy(m-ninc)=dsdy(m-ninc)+sn(m,kdir,2)*(s(m)+s(m-ninc))
                   dsdz(m-ninc)=dsdz(m-ninc)+sn(m,kdir,3)*(s(m)+s(m-ninc))
                enddo
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i1,j1  ,k)
             ind2 = indc(i1,j2m1,k)
!!$OMP SIMD
             do n=ind1,ind2,ncj
                m=n-n0c
                dsdx(m)=dsdx(m)-sn(m,kdir,1)*2*s(m-ninc)
                dsdy(m)=dsdy(m)-sn(m,kdir,2)*2*s(m-ninc)
                dsdz(m)=dsdz(m)-sn(m,kdir,3)*2*s(m-ninc)
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i2,j1  ,k)
             ind2 = indc(i2,j2m1,k)
!!$OMP SIMD
             do n=ind1,ind2,ncj
                m=n-n0c
                dsdx(m-ninc)=dsdx(m-ninc)+sn(m,kdir,1)*2*s(m)
                dsdy(m-ninc)=dsdy(m-ninc)+sn(m,kdir,2)*2*s(m)
                dsdz(m-ninc)=dsdz(m-ninc)+sn(m,kdir,3)*2*s(m)
             enddo
          enddo
!
       endif
!********************************************************************
    else   !schemas ponderes
!
!--------------------------------------------------------------------
!     direction k (a travers les facettes k=cste)
!--------------------------------------------------------------------
       if (equat(3:5).ne.'2dk') then
!
          kdir=3
          ninc=nck
!
          do k=k1p1,k2m1
             do j=j1,j2m1
                ind1 = indc(i1,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   sk1=(cmuk1(m)*s(m)+cmuk2(m)*s(m-ninc))*sn(m,kdir,1)
                   sk2=(cmuk1(m)*s(m)+cmuk2(m)*s(m-ninc))*sn(m,kdir,2)
                   sk3=(cmuk1(m)*s(m)+cmuk2(m)*s(m-ninc))*sn(m,kdir,3)
                   dsdx(m)=dsdx(m)-sk1
                   dsdy(m)=dsdy(m)-sk2
                   dsdz(m)=dsdz(m)-sk3
                enddo
                do n=ind1,ind2
                   m=n-n0c
                   sk1=(cmuk1(m)*s(m)+cmuk2(m)*s(m-ninc))*sn(m,kdir,1)
                   sk2=(cmuk1(m)*s(m)+cmuk2(m)*s(m-ninc))*sn(m,kdir,2)
                   sk3=(cmuk1(m)*s(m)+cmuk2(m)*s(m-ninc))*sn(m,kdir,3)
                   dsdx(m-ninc)=dsdx(m-ninc)+sk1
                   dsdy(m-ninc)=dsdy(m-ninc)+sk2
                   dsdz(m-ninc)=dsdz(m-ninc)+sk3
                enddo
             enddo
          enddo
!
          do j=j1,j2m1
             ind1 = indc(i1 ,j,k1)
             ind2 = indc(i2m1,j,k1)
             do n=ind1,ind2
                m=n-n0c
                dsdx(m)=dsdx(m)-sn(m,kdir,1)*2*s(m-ninc)
                dsdy(m)=dsdy(m)-sn(m,kdir,2)*2*s(m-ninc)
                dsdz(m)=dsdz(m)-sn(m,kdir,3)*2*s(m-ninc)
             enddo
          enddo
!
          do j=j1,j2m1
             ind1 = indc(i1 ,j,k2)
             ind2 = indc(i2m1,j,k2)
             do n=ind1,ind2
                m=n-n0c
                dsdx(m-ninc)=dsdx(m-ninc)+sn(m,kdir,1)*2*s(m)
                dsdy(m-ninc)=dsdy(m-ninc)+sn(m,kdir,2)*2*s(m)
                dsdz(m-ninc)=dsdz(m-ninc)+sn(m,kdir,3)*2*s(m)
             enddo
          enddo
!
       endif
!--------------------------------------------------------------
!     direction j (a travers les facettes j=cste)
!--------------------------------------------------------------
!
       if (equat(3:5).ne.'2dj') then
!
          kdir=2
          ninc=ncj
!
          do j=j1p1,j2m1
             do k=k1,k2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   sj1=(cmuj1(m)*s(m)+cmuj2(m)*s(m-ninc))*sn(m,kdir,1)
                   sj2=(cmuj1(m)*s(m)+cmuj2(m)*s(m-ninc))*sn(m,kdir,2)
                   sj3=(cmuj1(m)*s(m)+cmuj2(m)*s(m-ninc))*sn(m,kdir,3)
                   dsdx(m)=dsdx(m)-sj1
                   dsdy(m)=dsdy(m)-sj2
                   dsdz(m)=dsdz(m)-sj3
                enddo
                do n=ind1,ind2
                   m=n-n0c
                   sj1=(cmuj1(m)*s(m)+cmuj2(m)*s(m-ninc))*sn(m,kdir,1)
                   sj2=(cmuj1(m)*s(m)+cmuj2(m)*s(m-ninc))*sn(m,kdir,2)
                   sj3=(cmuj1(m)*s(m)+cmuj2(m)*s(m-ninc))*sn(m,kdir,3)
                   dsdx(m-ninc)=dsdx(m-ninc)+sj1
                   dsdy(m-ninc)=dsdy(m-ninc)+sj2
                   dsdz(m-ninc)=dsdz(m-ninc)+sj3
                enddo
             enddo
          enddo
!
          do k=k1,k2m1
             ind1=indc(i1  ,j1,k)
             ind2=indc(i2m1,j1,k)
             do n=ind1,ind2
                m=n-n0c
                dsdx(m)=dsdx(m)-sn(m,kdir,1)*2*s(m-ninc)
                dsdy(m)=dsdy(m)-sn(m,kdir,2)*2*s(m-ninc)
                dsdz(m)=dsdz(m)-sn(m,kdir,3)*2*s(m-ninc)
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i1  ,j2,k)
             ind2 = indc(i2m1,j2,k)
             do n=ind1,ind2
                m=n-n0c
                dsdx(m-ninc)=dsdx(m-ninc)+sn(m,kdir,1)*2*s(m)
                dsdy(m-ninc)=dsdy(m-ninc)+sn(m,kdir,2)*2*s(m)
                dsdz(m-ninc)=dsdz(m-ninc)+sn(m,kdir,3)*2*s(m)
             enddo
          enddo
!
       endif
!---------------------------------------------------------------------
!     direction i (a travers les facettes i=cste)
!---------------------------------------------------------------------
       if (equat(3:5).ne.'2di') then
!
          kdir=1
          ninc=nci
!
          do k=k1,k2m1
             do j=j1,j2m1
                ind1 = indc(i1p1,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   si1=(cmui1(m)*s(m)+cmui2(m)*s(m-ninc))*sn(m,kdir,1)
                   si2=(cmui1(m)*s(m)+cmui2(m)*s(m-ninc))*sn(m,kdir,2)
                   si3=(cmui1(m)*s(m)+cmui2(m)*s(m-ninc))*sn(m,kdir,3)
                   dsdx(m)=dsdx(m)-si1
                   dsdy(m)=dsdy(m)-si2
                   dsdz(m)=dsdz(m)-si3
                enddo
                do n=ind1,ind2
                   m=n-n0c
                   si1=(cmui1(m)*s(m)+cmui2(m)*s(m-ninc))*sn(m,kdir,1)
                   si2=(cmui1(m)*s(m)+cmui2(m)*s(m-ninc))*sn(m,kdir,2)
                   si3=(cmui1(m)*s(m)+cmui2(m)*s(m-ninc))*sn(m,kdir,3)
                   dsdx(m-ninc)=dsdx(m-ninc)+si1
                   dsdy(m-ninc)=dsdy(m-ninc)+si2
                   dsdz(m-ninc)=dsdz(m-ninc)+si3
                enddo
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i1,j1  ,k)
             ind2 = indc(i1,j2m1,k)
!!$OMP SIMD
             do n=ind1,ind2,ncj
                m=n-n0c
                dsdx(m)=dsdx(m)-sn(m,kdir,1)*2*s(m-ninc)
                dsdy(m)=dsdy(m)-sn(m,kdir,2)*2*s(m-ninc)
                dsdz(m)=dsdz(m)-sn(m,kdir,3)*2*s(m-ninc)
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i2,j1  ,k)
             ind2 = indc(i2,j2m1,k)
!!$OMP SIMD
             do n=ind1,ind2,ncj
                m=n-n0c
                dsdx(m-ninc)=dsdx(m-ninc)+sn(m,kdir,1)*2*s(m)
                dsdy(m-ninc)=dsdy(m-ninc)+sn(m,kdir,2)*2*s(m)
                dsdz(m-ninc)=dsdz(m-ninc)+sn(m,kdir,3)*2*s(m)
             enddo
          enddo
!
       endif
!
    endif  ! fin test sur schema
!-----------------------------------------------------------------
!-     calcul du gradient :
!----------------------------------------------------------------
!
    ind1 = indc(i1  ,j1  ,k1  )
    ind2 = indc(i2m1,j2m1,k2m1)
    do n=ind1,ind2
       m=n-n0c
!       le coefficient 1/2 provient de la moyenne de vx,vy,vz ou t
       ts=sign(0.5,-vol(n))
       vols=(0.5+ts)*eps+(0.5-ts)*vol(n)
       c0=0.5/vols
!        c0=0.5/vol(n)
       dsdx(m)=dsdx(m)*c0
       dsdy(m)=dsdy(m)*c0
       dsdz(m)=dsdz(m)*c0
    enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
  integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
    function    inc(id,jd,kd)
      implicit none
  integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine teq_grads
end module mod_teq_grads
