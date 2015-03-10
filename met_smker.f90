module mod_met_smker
implicit none
contains
      subroutine met_smker( &
                 l, &
                 s,cfke, &
                 tprod,tsv6,tsv7, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz)
!
!***********************************************************************
!
!_DA  decembre 2000 - Eric Goncalves
!
!     calcul du terme source pour modele k-eps avec correction de Durbin
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
implicit none
integer :: inc
integer :: indc
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: s
double precision :: cfke
double precision :: tprod
double precision :: tsv6
double precision :: tsv7
double precision :: dvxx
double precision :: dvxy
double precision :: dvxz
double precision :: dvyx
double precision :: dvyy
double precision :: dvyz
double precision :: dvzx
double precision :: dvzy
double precision :: dvzz
double precision :: as
double precision :: b
double precision :: c1f1
double precision :: c2f2
double precision :: d
double precision :: esk
double precision :: g
integer :: i1
integer :: i1m1
integer :: i1p1
integer :: i2
integer :: i2m1
integer :: j1
integer :: j1p1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k1p1
integer :: k2
integer :: k2m1
integer :: m
integer :: n
integer :: n0c
integer :: nci
integer :: nid
integer :: nijd
integer :: njd
double precision :: rdelta
double precision :: roe2sk
double precision :: x
double precision :: xl1
double precision :: xl2
!
!-----------------------------------------------------------------------
!
      dimension s(ip11,ip60)
      dimension tsv6(ip12),tsv7(ip12)
      dimension tprod(ip00)
      dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
                dvyx(ip00),dvyy(ip00),dvyz(ip00), &
                dvzx(ip00),dvzy(ip00),dvzz(ip00)
      dimension cfke(ip13)
!
      indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
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
      i1m1=i1-1
!
      nci  = inc(1,0,0)
!
      as=1./sqrt(cmu)
      do k=k1,k2m1
        do j=j1,j2m1
          n=indc(i1m1,j,k)
          do i=i1,i2m1
            n=n+nci
            m=n-n0c
            esk=s(n,7)/s(n,6)
            roe2sk=(s(n,7)*s(n,7))/s(n,6)
!            ss=s(n,6)*sqrt(4.*(dvxx(m)**2+dvyy(m)**2+
!     &            dvzz(m)**2)/3. + (dvzy(m)+dvyz(m))**2
!     &              + (dvxz(m)+dvzx(m))**2
!     &              + (dvyx(m)+dvxy(m))**2)/s(n,7)
!            alpha=min(1.,1./(as*cmu*ss))
            tsv6(n)=tprod(m)-s(n,7)
            tsv7(n)=cke1*tprod(m)*esk-cke2*roe2sk
!           tsv7(n)=cke1*tprod(m)*esk/alpha -cke2*roe2sk
!           rayon spectral
            c1f1=cke1
            c2f2=cke2
            d=0.
            x=c1f1*tprod(m)-c2f2*s(n,7)
            g=max(0.,d*(d-2.*(x-c2f2*s(n,7))))
            rdelta=sqrt((x-c2f2*s(n,7))**2+4*x*s(n,7)+g)
            b=x-c2f2*s(n,7)-d
            xl1=0.5*(b+rdelta)/s(n,6)
            xl2=0.5*(b-rdelta)/s(n,6)
            cfke(n)=max(abs(xl1),abs(xl2))
          enddo
        enddo
      enddo
!
      return
      end
end module
