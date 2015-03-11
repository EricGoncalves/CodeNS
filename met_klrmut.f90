module mod_met_klrmut
implicit none
contains
      subroutine met_klrmut( &
                 l, &
                 sn,vol,t, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 dist,v,mu,mut, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DA  DATE_C :  janvier 2010 - Jean Decaix
!
!     ACT
!_A   Calcul de la viscosite turbulente a partir de k et l
!_A   Modele de Smith - modele avec correction de realisabilite
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
double precision :: dist
double precision :: v
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
double precision :: a1
double precision :: c1
double precision :: c14
double precision :: c2
double precision :: c22
double precision :: coef1
double precision :: coef2
double precision :: exp2x
double precision :: f1
double precision :: f2
double precision :: fmu
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
integer :: ncj
integer :: nck
integer :: nid
integer :: nijd
integer :: njd
double precision :: r2sb1
double precision :: rack
double precision :: ss2
double precision :: tdef
double precision :: xl2
double precision :: xxi
double precision :: xxi2
double precision :: zeta
!
!-----------------------------------------------------------------------
!
      real mu,mut
      dimension mu(ip12),mut(ip12),dist(ip12)
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
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd= nid*njd
!
      nci = inc(1,0,0)
      ncj = inc(0,1,0)
      nck = inc(0,0,1)
!
!     notation pour les constantes
!
!     B1 <-> cklb1
!     E2 <-> ckle2
!     sigma_k <-> sigmak
!     sigma_l <-> sigmal
!
      c1=25.5
      c2=2.
!
      r2sb1 =sqrt(2.)/cklb1**(1./3.)
      c14=c1**4
      c22=c2**2
!
!     correction SST
!      a1=sqrt(cmukl)      !classique a1=0.3
      a1=0.3   

!
!     Calcul du gradient de la vitesse. Les tableaux ont ete reutilises
!     dans la phase implicite.
!
        call teq_gradv( &
             l, &
             sn, &
             vol,v, &
             t , &
             dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
             cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
        do k=k1,k2m1
         do j=j1,j2m1
          n=ind(i1-1,j,k)
          do i=i1,i2m1
           n=n+nci
           m=n-n0
           rack=sqrt(v(n,6)/v(n,1))
           xxi=r2sb1*rack*v(n,7)/mu(n)
           xxi2=xxi*xxi
           xl2=(v(n,7)/(xkappa*v(n,1)*dist(n)))**2
           f1=exp(-50.*xl2)
           ss2=(c22+xxi2)*xxi2
           fmu=((c14*f1+ss2)/(c14+ss2))**0.25
!
           coef1=2.*r2sb1*v(n,7)/(cmukl*dist(n)*v(n,1))
           coef2=500.*r2sb1*mu(n)*v(n,7)/(rack*v(n,1)**2*dist(n)**2)
           tdef = sqrt( (dvzy(m)+dvyz(m))**2 &
                       +(dvxz(m)+dvzx(m))**2 &
                       +(dvyx(m)+dvxy(m))**2)
           zeta=max(coef1,coef2)
           exp2x=exp(min(2.*zeta**2,25.))
           f2=(exp2x-1.)/(exp2x+1.)
           mut(n)=v(n,6)/max(rack*v(n,1)/(fmu*v(n,7)*r2sb1),tdef*f2/a1)
          enddo
         enddo
        enddo
!
      return
      end subroutine
end module
