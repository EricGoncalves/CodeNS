module mod_met_laplaciens
  implicit none
contains
  subroutine met_laplaciens( &
       l, &
       equat, &
       sn,lgsnlt, &
       vol, &
       v, &
       dsx,dsy,dsz, &
       dsd2x,dsd2y,dsd2z, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DATE: Mars 2010 - Jean DECAIX / LEGI
!
!     ACT     calcul du laplacien de la vitesse
!
!       dsx      (ip 00)       composante suivant x du gradient de la vitesse
!       dsy      (ip 00)       composante suivant y du gradient de la vitesse
!       dsz      (ip 00)       composante suivant z du gradient de la vitesse
!       dsd2x    (ip 00)       dérivée de ds/x
!       dsd2y    (ip 00)       dérivéé de ds/y
!       dsd2z    (ip 00)       dérivée de ds/z
!       v        (ip 11 ip 60) variables à l'instant n
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
    double precision ::    c0,cmui1,cmui2,cmuj1,cmuj2
    double precision :: cmuk1,cmuk2,dsd2x,dsd2y,dsd2z
    double precision ::   dsx,  dsy,  dsz,  eps,  si1
    double precision ::   si2,  si3,  sj1,  sj2,  sj3
    double precision ::   sk1,  sk2,  sk3,   sn,   ts
    double precision ::     v,  vol, vols
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
    dimension sn(lgsnlt,nind,ndir), &
         vol(ip11)
    dimension dsx(ip00),dsy(ip00),dsz(ip00), &
              dsd2x(ip00),dsd2y(ip00),dsd2z(ip00)
    dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
              cmuk1(ip21),cmuk2(ip21)
    dimension v(ip11,ip60)

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
       dsd2x(m)=0.
       dsd2y(m)=0.
       dsd2z(m)=0.
    enddo
!
!cc   test sur le schema numerique
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
!!$OMP SIMD
                do n=ind1,ind2
                   m=n-n0c
                   dsd2x(m)=dsd2x(m)-sn(m,kdir,1)*(dsx(m)+dsx(m-ninc))
                   dsd2y(m)=dsd2y(m)-sn(m,kdir,2)*(dsy(m)+dsy(m-ninc))
                   dsd2z(m)=dsd2z(m)-sn(m,kdir,3)*(dsz(m)+dsz(m-ninc))
                   dsd2x(m-ninc)=dsd2x(m-ninc)+sn(m,kdir,1)*(dsx(m)+dsx(m-ninc))
                   dsd2y(m-ninc)=dsd2y(m-ninc)+sn(m,kdir,2)*(dsy(m)+dsy(m-ninc))
                   dsd2z(m-ninc)=dsd2z(m-ninc)+sn(m,kdir,3)*(dsz(m)+dsz(m-ninc))
                enddo
             enddo
          enddo
!
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k1)
             ind2 = indc(i2m1,j,k1)
!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                dsd2x(m)=dsd2x(m)-sn(m,kdir,1)*2*dsx(m-ninc)
                dsd2y(m)=dsd2y(m)-sn(m,kdir,2)*2*dsy(m-ninc)
                dsd2z(m)=dsd2z(m)-sn(m,kdir,3)*2*dsz(m-ninc)
             enddo
          enddo
!
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k2)
             ind2 = indc(i2m1,j,k2)
!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                dsd2x(m-ninc)=dsd2x(m-ninc)+sn(m,kdir,1)*2*dsx(m)
                dsd2y(m-ninc)=dsd2y(m-ninc)+sn(m,kdir,2)*2*dsy(m)
                dsd2z(m-ninc)=dsd2z(m-ninc)+sn(m,kdir,3)*2*dsz(m)
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
!!$OMP SIMD
                do n=ind1,ind2
                   m=n-n0c
                   dsd2x(m)=dsd2x(m)-sn(m,kdir,1)*(dsx(m)+dsy(m-ninc))
                   dsd2y(m)=dsd2y(m)-sn(m,kdir,2)*(dsy(m)+dsy(m-ninc))
                   dsd2z(m)=dsd2z(m)-sn(m,kdir,3)*(dsz(m)+dsz(m-ninc))
                   dsd2x(m-ninc)=dsd2x(m-ninc)+sn(m,kdir,1)*(dsx(m)+dsx(m-ninc))
                   dsd2y(m-ninc)=dsd2y(m-ninc)+sn(m,kdir,2)*(dsy(m)+dsy(m-ninc))
                   dsd2z(m-ninc)=dsd2z(m-ninc)+sn(m,kdir,3)*(dsz(m)+dsz(m-ninc))
                enddo
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i1  ,j1,k)
             ind2 = indc(i2m1,j1,k)
!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                dsd2x(m)=dsd2x(m)-sn(m,kdir,1)*2*dsx(m-ninc)
                dsd2y(m)=dsd2y(m)-sn(m,kdir,2)*2*dsy(m-ninc)
                dsd2z(m)=dsd2z(m)-sn(m,kdir,3)*2*dsz(m-ninc)
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i1  ,j2,k)
             ind2 = indc(i2m1,j2,k)
!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                dsd2x(m-ninc)=dsd2x(m-ninc)+sn(m,kdir,1)*2*dsx(m)
                dsd2y(m-ninc)=dsd2y(m-ninc)+sn(m,kdir,2)*2*dsy(m)
                dsd2z(m-ninc)=dsd2z(m-ninc)+sn(m,kdir,3)*2*dsz(m)
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
!!$OMP SIMD
                do n=ind1,ind2
                   m=n-n0c
                   dsd2x(m)=dsd2x(m)-sn(m,kdir,1)*(dsx(m)+dsx(m-ninc))
                   dsd2y(m)=dsd2y(m)-sn(m,kdir,2)*(dsy(m)+dsy(m-ninc))
                   dsd2z(m)=dsd2z(m)-sn(m,kdir,3)*(dsz(m)+dsz(m-ninc))
                   dsd2x(m-ninc)=dsd2x(m-ninc)+sn(m,kdir,1)*(dsx(m)+dsx(m-ninc))
                   dsd2y(m-ninc)=dsd2y(m-ninc)+sn(m,kdir,2)*(dsy(m)+dsy(m-ninc))
                   dsd2z(m-ninc)=dsd2z(m-ninc)+sn(m,kdir,3)*(dsz(m)+dsz(m-ninc))
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
                dsd2x(m)=dsd2x(m)-sn(m,kdir,1)*2*dsx(m-ninc)
                dsd2y(m)=dsd2y(m)-sn(m,kdir,2)*2*dsy(m-ninc)
                dsd2z(m)=dsd2z(m)-sn(m,kdir,3)*2*dsz(m-ninc)
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i2,j1  ,k)
             ind2 = indc(i2,j2m1,k)
!!$OMP SIMD
             do n=ind1,ind2,ncj
                m=n-n0c
                dsd2x(m-ninc)=dsd2x(m-ninc)+sn(m,kdir,1)*2*dsx(m)
                dsd2y(m-ninc)=dsd2y(m-ninc)+sn(m,kdir,2)*2*dsy(m)
                dsd2z(m-ninc)=dsd2z(m-ninc)+sn(m,kdir,3)*2*dsz(m)
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
!!$OMP SIMD
                do n=ind1,ind2
                   m=n-n0c
                   sk1=(cmuk1(n)*dsx(m)+cmuk2(n)*dsx(m-ninc))*sn(m,kdir,1)
                   sk2=(cmuk1(n)*dsy(m)+cmuk2(n)*dsy(m-ninc))*sn(m,kdir,2)
                   sk3=(cmuk1(n)*dsz(m)+cmuk2(n)*dsz(m-ninc))*sn(m,kdir,3)
                   dsd2x(m)=dsd2x(m)-sk1
                   dsd2y(m)=dsd2y(m)-sk2
                   dsd2z(m)=dsd2z(m)-sk3
                   dsd2x(m-ninc)=dsd2x(m-ninc)+sk1
                   dsd2y(m-ninc)=dsd2y(m-ninc)+sk2
                   dsd2z(m-ninc)=dsd2z(m-ninc)+sk3
                enddo
             enddo
          enddo
!
          do j=j1,j2m1
             ind1 = indc(i1 ,j,k1)
             ind2 = indc(i2m1,j,k1)
!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                dsd2x(m)=dsd2x(m)-sn(m,kdir,1)*2*dsx(m-ninc)
                dsd2y(m)=dsd2y(m)-sn(m,kdir,2)*2*dsy(m-ninc)
                dsd2z(m)=dsd2z(m)-sn(m,kdir,3)*2*dsz(m-ninc)
             enddo
          enddo
!
          do j=j1,j2m1
             ind1 = indc(i1 ,j,k2)
             ind2 = indc(i2m1,j,k2)
!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                dsd2x(m-ninc)=dsd2x(m-ninc)+sn(m,kdir,1)*2*dsx(m)
                dsd2y(m-ninc)=dsd2y(m-ninc)+sn(m,kdir,2)*2*dsy(m)
                dsd2z(m-ninc)=dsd2z(m-ninc)+sn(m,kdir,3)*2*dsz(m)
             enddo
          enddo
!
       endif
!--------------------------------------------------------------
!     direction j (a travers les facettes j=cste)
!--------------------------------------------------------------
       if (equat(3:5).ne.'2dj') then
!
          kdir=2
          ninc=ncj
!
          do j=j1p1,j2m1
             do k=k1,k2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
!!$OMP SIMD
                do n=ind1,ind2
                   m=n-n0c
                   sj1=(cmuj1(n)*dsx(m)+cmuj2(n)*dsx(m-ninc))*sn(m,kdir,1)
                   sj2=(cmuj1(n)*dsy(m)+cmuj2(n)*dsy(m-ninc))*sn(m,kdir,2)
                   sj3=(cmuj1(n)*dsz(m)+cmuj2(n)*dsz(m-ninc))*sn(m,kdir,3)
                   dsd2x(m)=dsd2x(m)-sj1
                   dsd2y(m)=dsd2y(m)-sj2
                   dsd2z(m)=dsd2z(m)-sj3
                   dsd2x(m-ninc)=dsd2x(m-ninc)+sj1
                   dsd2y(m-ninc)=dsd2y(m-ninc)+sj2
                   dsd2z(m-ninc)=dsd2z(m-ninc)+sj3
                enddo
             enddo
          enddo
!
          do k=k1,k2m1
             ind1=indc(i1  ,j1,k)
             ind2=indc(i2m1,j1,k)
!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                dsd2x(m)=dsd2x(m)-sn(m,kdir,1)*2*dsx(m-ninc)
                dsd2y(m)=dsd2y(m)-sn(m,kdir,2)*2*dsy(m-ninc)
                dsd2z(m)=dsd2z(m)-sn(m,kdir,3)*2*dsz(m-ninc)
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i1  ,j2,k)
             ind2 = indc(i2m1,j2,k)
!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0c
                dsd2x(m-ninc)=dsd2x(m-ninc)+sn(m,kdir,1)*2*dsx(m)
                dsd2y(m-ninc)=dsd2y(m-ninc)+sn(m,kdir,2)*2*dsy(m)
                dsd2z(m-ninc)=dsd2z(m-ninc)+sn(m,kdir,3)*2*dsz(m)
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
!!$OMP SIMD
                do n=ind1,ind2
                   m=n-n0c
                   si1=(cmui1(n)*dsx(m)+cmui2(n)*dsx(m-ninc))*sn(m,kdir,1)
                   si2=(cmui1(n)*dsy(m)+cmui2(n)*dsy(m-ninc))*sn(m,kdir,2)
                   si3=(cmui1(n)*dsz(m)+cmui2(n)*dsz(m-ninc))*sn(m,kdir,3)
                   dsd2x(m)=dsd2x(m)-si1
                   dsd2y(m)=dsd2y(m)-si2
                   dsd2z(m)=dsd2z(m)-si3
                   dsd2x(m-ninc)=dsd2x(m-ninc)+si1
                   dsd2y(m-ninc)=dsd2y(m-ninc)+si2
                   dsd2z(m-ninc)=dsd2z(m-ninc)+si3
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
                dsd2x(m)=dsd2x(m)-sn(m,kdir,1)*2*dsx(m-ninc)
                dsd2y(m)=dsd2y(m)-sn(m,kdir,2)*2*dsy(m-ninc)
                dsd2z(m)=dsd2z(m)-sn(m,kdir,3)*2*dsz(m-ninc)
             enddo
          enddo
!
          do k=k1,k2m1
             ind1 = indc(i2,j1  ,k)
             ind2 = indc(i2,j2m1,k)
!!$OMP SIMD
             do n=ind1,ind2,ncj
                m=n-n0c
                dsd2x(m-ninc)=dsd2x(m-ninc)+sn(m,kdir,1)*2*dsx(m)
                dsd2y(m-ninc)=dsd2y(m-ninc)+sn(m,kdir,2)*2*dsy(m)
                dsd2z(m-ninc)=dsd2z(m-ninc)+sn(m,kdir,3)*2*dsz(m)
             enddo
          enddo
!
       endif
!
    endif  ! fin test sur schema
!-----------------------------------------------------------------
!-     calcul du laplacien :
!----------------------------------------------------------------
!
    ind1 = indc(i1  ,j1  ,k1  )
    ind2 = indc(i2m1,j2m1,k2m1)
    do n=ind1,ind2
       m=n-n0c
!       le coefficient 1/2 provient de la moyenne de vx,vy,vz ou t
       ts=sign(0.5D0,-vol(n))
       vols=(0.5+ts)*eps+(0.5-ts)*vol(n)
       c0=0.5/vols
!        c0=0.5/vol(n)
       dsd2x(m)=dsd2x(m)*c0
       dsd2y(m)=dsd2y(m)*c0
       dsd2z(m)=dsd2z(m)*c0
!      laplacien de ds
       dsd2x(m)=dsd2x(m)+dsd2y(m)+dsd2z(m)
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
  end subroutine met_laplaciens
end module mod_met_laplaciens
