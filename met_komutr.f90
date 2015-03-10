module mod_met_komutr
implicit none
contains
      subroutine met_komutr( &
                 l, &
                 sn,vol,t, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 dist,v,mu,mut, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_H   DATE_C : septembre 2002  - AUTEUR: Eric Goncalves - LEGI
!_H
!     ACT
!_A   Calcul de la viscosite turbulente a partir de k et omega
!_A   Modele de Menter avec conditions de realisabilite de Menter
!_A   avec ou sans correction SST.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!
!     OUT
!_O    mut        : arg real(ip12      ) ; viscosite turbulente
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
use mod_teq_gradv
implicit none
integer :: ind
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
double precision :: alpha
double precision :: betae2
double precision :: coef1
double precision :: coef2
double precision :: exp2x
double precision :: f2
integer :: i1
integer :: i2
integer :: i2m1
integer :: ind1
integer :: ind2
integer :: j1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k2
integer :: k2m1
integer :: m
integer :: n
integer :: n0
integer :: nid
integer :: nijd
integer :: njd
double precision :: omeg
double precision :: rcmu
double precision :: rmut
double precision :: rota
double precision :: smut
double precision :: ss
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
      dimension sn(ip31*ndir), &
                vol(ip11)
      dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
      ind(i,j,k)=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
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
      rcmu=1./sqrt(0.09)         !correspond a c=0.52
      betae2=betae/2.
      a1=sqrt(betae)
!
!       Calcul du gradient de la vitesse
        call teq_gradv( &
             l, &
             sn, &
             vol,v, &
             t , &
             dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
             cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
      if(komsst.eq.0) then
!       Modele sans correction SST
!
        do k=k1,k2m1
         do j=j1,j2m1
          ind1 = ind(i1  ,j,k)
          ind2 = ind(i2m1,j,k)
          do n=ind1,ind2
           m=n-n0
           ss=v(n,1)*sqrt(4.*(dvxx(m)**2+dvyy(m)**2+ &
              dvzz(m)**2)/3. + (dvzy(m)+dvyz(m))**2 &
                             + (dvxz(m)+dvzx(m))**2 &
               + (dvyx(m)+dvxy(m))**2)/v(n,7)
           alpha=min(1.,rcmu/ss)
           mut(n)=v(n,1)*alpha*v(n,6)/v(n,7)
          enddo
         enddo
        enddo
!
      else if(komsst.eq.1) then
!       Modele avec correction SST
!
        do k=k1,k2m1
         do j=j1,j2m1
          ind1 = ind(i1  ,j,k)
          ind2 = ind(i2m1,j,k)
          do n=ind1,ind2
           m=n-n0
           omeg=v(n,7)/v(n,1)
           coef1=sqrt(v(n,6)/v(n,1))/(betae2*omeg*dist(n))
           coef2=500.*mu(n)/(v(n,7)*dist(n)**2)
           rota = sqrt( (dvzy(m)-dvyz(m))**2 &
                +(dvxz(m)-dvzx(m))**2 &
                +(dvyx(m)-dvxy(m))**2)
           zeta=max(coef1,coef2)
           exp2x=exp(min(2.*zeta**2,25.))
           f2=(exp2x-1.)/(exp2x+1)
           smut=v(n,6)/max(omeg,rota*f2/a1)
           ss=v(n,1)*sqrt(4.*(dvxx(m)**2+dvyy(m)**2+ &
              dvzz(m)**2)/3. + (dvzy(m)+dvyz(m))**2 &
                             + (dvxz(m)+dvzx(m))**2 &
               + (dvyx(m)+dvxy(m))**2)/v(n,7)
           alpha=min(1.,rcmu/ss)
           rmut=v(n,1)*alpha*v(n,6)/v(n,7)
           mut(n)=min(smut,rmut)
          enddo
         enddo
        enddo
      endif
!
      return
      end
end module
