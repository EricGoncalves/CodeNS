module mod_teq_gradv
implicit none
contains
      subroutine teq_gradv( &
                 l, &
                 sn, &
                 vol,t, &
                 s, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT  calcul du tenseur gradient de la vitesse
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use chainecarac
use mod_teq_grads
implicit none
integer :: indc
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: sn
double precision :: vol
double precision :: t
double precision :: s
double precision :: dvxx
double precision :: dvxy
double precision :: dvxz
double precision :: dvyx
double precision :: dvyy
double precision :: dvyz
double precision :: dvzx
double precision :: dvzy
double precision :: dvzz
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
integer :: i1
integer :: i1m1
integer :: i2
integer :: i2m1
integer :: imax
integer :: imin
integer :: ind1
integer :: ind2
integer :: j1
integer :: j1m1
integer :: j2
integer :: j2m1
integer :: jmax
integer :: jmin
integer :: k1
integer :: k1m1
integer :: k2
integer :: k2m1
integer :: kmax
integer :: kmin
integer :: lgsnlt
integer :: m
integer :: n
integer :: n0
integer :: nid
integer :: nijd
integer :: njd
integer :: npsn
!
!-----------------------------------------------------------------------
!
      dimension t(ip11,ip60),vol(ip11)
      dimension sn(ip31*ndir)
      dimension s(ip00),dvxx(ip00),dvxy(ip00),dvxz(ip00), &
                dvyx(ip00),dvyy(ip00),dvyz(ip00), &
                dvzx(ip00),dvzy(ip00),dvzz(ip00)
      dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
      indc(i,j,k)=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
!
      n0=npc(l)
      i1=ii1(l)
      i2=ii2(l)
      j1=jj1(l)
      j2=jj2(l)
      k1=kk1(l)
      k2=kk2(l)
!
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
!
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd= nid*njd
!
      imin=i1m1
      imax=i2
      jmin=j1m1
      jmax=j2
      kmin=k1m1
      kmax=k2
!
      if (equat(3:5).eq.'2di') then
         imin=i1
         imax=i2m1
      endif
      if (equat(3:5).eq.'2dj') then
         jmin=j1
         jmax=j2m1
      endif
      if (equat(3:5).eq.'2dk') then
         kmin=k1
         kmax=k2m1
      endif
!
!-----calcul de grad(vx) a l'instant n
!
      do k=kmin,kmax
       do j=jmin,jmax
        ind1=indc(imin,j,k)
        ind2=indc(imax,j,k)
        do n=ind1,ind2
         t(n,1) = max(t(n,1),1.e-20)
        enddo
       enddo
      enddo
!
      do k=kmin,kmax
       do j=jmin,jmax
        ind1=indc(imin,j,k)
        ind2=indc(imax,j,k)
        do n=ind1,ind2
          m=n-n0
          s(m)= t(n,2)/t(n,1)
        enddo
       enddo
      enddo
!
      npsn  =ndir*npfb(l)+1
      lgsnlt=nnn(l)
!
      call teq_grads( &
                 l, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 s, &
                 dvxx,dvxy,dvxz, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!-----calcul de grad(vy) a l'instant n
!
      do k=kmin,kmax
        do j=jmin,jmax
          ind1=indc(imin,j,k)
          ind2=indc(imax,j,k)
          do n=ind1,ind2
            m=n-n0
            s(m)= t(n,3)/t(n,1)
          enddo
        enddo
      enddo
!
      npsn  =ndir*npfb(l)+1
      lgsnlt=nnn(l)
!
      call teq_grads( &
                 l, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 s, &
                 dvyx,dvyy,dvyz, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!-----calcul de grad(vz) a l'instant n
!
      do k=kmin,kmax
        do j=jmin,jmax
          ind1=indc(imin,j,k)
          ind2=indc(imax,j,k)
          do n=ind1,ind2
            m=n-n0
            s(m)= t(n,4)/t(n,1)
          enddo
        enddo
      enddo
!
      npsn  =ndir*npfb(l)+1
      lgsnlt=nnn(l)
!
      call teq_grads( &
                 l, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 s, &
                 dvzx,dvzy,dvzz, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
      return
      end subroutine
end module
