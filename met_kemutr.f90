module mod_met_kemutr
implicit none
contains
      subroutine met_kemutr( &
                 l,ncyc, &
                 sn,vol,t, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 v,mu,mut, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DA : decembre 2000 - Eric Goncalves
!
!     ACT
!      Calcul de la viscosite turbulente mut avec correction de 
!       realisabilite de Durbin.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
use mod_teq_gradv
implicit none
integer :: inc
integer :: ind
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: l
integer :: ncyc
double precision :: sn
double precision :: vol
double precision :: t
double precision :: dvxx
double precision :: dvxy
double precision :: dvxz
double precision :: dvyx
double precision :: dvyy
double precision :: dvyz
double precision :: dvzx
double precision :: dvzy
double precision :: dvzz
double precision :: v
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
integer :: i1
integer :: i2
integer :: i2m1
integer :: j1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k2
integer :: k2m1
integer :: m
integer :: n
integer :: n0
integer :: nci
integer :: nid
integer :: nijd
integer :: njd
!
!-----------------------------------------------------------------------
!
      double precision mu,mut,cmuv,as,ss,retur
      dimension mu(ip12),mut(ip12)
      dimension v(ip11,ip60)
      dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
                dvyx(ip00),dvyy(ip00),dvyz(ip00), &
                dvzx(ip00),dvzy(ip00),dvzz(ip00),t(ip00)
      dimension sn(ip31*ndir),vol(ip11)
      dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
      ind(i,j,k)=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
!
      n0=npc(l)
      i1=ii1(l)
      i2=ii2(l)
      j1=jj1(l)
      j2=jj2(l)
      k1=kk1(l)
      k2=kk2(l)
!
!
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd= nid*njd
!
      nci = inc(1,0,0)
!
!       Calcul du gradient de la vitesse. Les tableaux ont ete reutilises
!       dans la phase implicite.
        call teq_gradv( &
             l, &
             sn, &
             vol,v, &
             t, &
             dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
             cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!        as=sqrt(cmu)         !correspond a c=0.52
        as=0.3/sqrt(3.)      !correspond a c=0.3
        do k=k1,k2m1
          do j=j1,j2m1
            n=ind(i1-1,j,k)
            do i=i1,i2m1
              n=n+nci
              m=n-n0
              retur=(v(n,6)**2)/(v(n,7)*mu(n))
              ss=v(n,6)*sqrt(4.*(dvxx(m)**2+dvyy(m)**2+ dvzz(m)**2)/3.  &
          + (dvzy(m)+dvyz(m))**2 &
                   + (dvxz(m)+dvzx(m))**2 &
                   + (dvyx(m)+dvxy(m))**2)/v(n,7)
              cmuv=min(cmu,as/ss)
              mut(n)=cmuv*retur*mu(n)
            enddo
          enddo
        enddo
!
      return
      end subroutine
end module
