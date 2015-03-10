module mod_met_brkl
implicit none
contains
      subroutine met_brkl( &
                 l, &
                 mu,mut,v, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz)
!
!***********************************************************************
!
!     ACT
!_A   initialisation des grandeurs k et l pour le modele de Smith
!_A   a partir de mu_t (u'v'/k=0.3)
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
implicit none
integer :: indc
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: v
double precision :: dvxx
double precision :: dvxy
double precision :: dvxz
double precision :: dvyx
double precision :: dvyy
double precision :: dvyz
double precision :: dvzx
double precision :: dvzy
double precision :: dvzz
integer :: i1
integer :: i1m1
integer :: i2
integer :: i2m1
integer :: ind1
integer :: ind2
integer :: j1
integer :: j1m1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k1m1
integer :: k2
integer :: k2m1
integer :: m
integer :: n
integer :: n0
integer :: nid
integer :: nijd
integer :: njd
!
!-----------------------------------------------------------------------
!
      real mu,mut,mutmx,omeg,b1s2,cmu1
      dimension v(ip11,ip60)
      dimension mu(ip12),mut(ip12)
      dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
                dvyx(ip00),dvyy(ip00),dvyz(ip00), &
                dvzx(ip00),dvzy(ip00),dvzz(ip00)
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
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd= nid*njd
!
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
!     notation pour les constantes
!
!     B1 <-> cklb1
!     E2 <-> ckle2
!     sigma_k <-> sigmak
!     sigma_l <-> sigmal
!
      b1s2=cklb1*(1./3.)/sqrt(2.)
      cmu1=1./sqrt(cmukl)
!
      do k=k1,k2m1
       do j=j1,j2m1
        ind1=indc(i1,j,k)
        ind2=indc(i2m1,j,k)
        do n=ind1,ind2
         m=n-n0
         omeg = sqrt((dvzy(m)-dvyz(m))**2 &
                    +(dvxz(m)-dvzx(m))**2 &
                    +(dvyx(m)-dvxy(m))**2)
         mutmx=max(mut(n),1.e-3*mu(n))
         v(n,6)=max(cmu1*mutmx*omeg,epsk)
         v(n,7)=b1s2*mutmx/sqrt(v(n,6)/v(n,1))
        enddo
       enddo
      enddo
!
      return
      end
end module
