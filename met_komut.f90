module mod_met_komut
use mod_teq_gradv
implicit none
contains
      subroutine met_komut( &
                 l, &
                 sn,vol,t, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 dist,v,mu,mut, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT
!_A   Calcul de la viscosite turbulente a partir de k et omega
!_A   Modele de Wilcox et Menter
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!
!     OUT
!_O    mut        : arg real(ip12      ) ; viscosite turbulente
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
      use proprieteflu
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
double precision :: betae2
double precision :: coef1
double precision :: coef2
double precision :: exp2x
double precision :: f2
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
double precision :: omeg
double precision :: rota
double precision :: zeta
!
!-----------------------------------------------------------------------
!
      double precision mu,mut
      dimension mu(ip12),mut(ip12),dist(ip12),vol(ip11)
      dimension v(ip11,ip60)
      dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
                dvyx(ip00),dvyy(ip00),dvyz(ip00), &
                dvzx(ip00),dvzy(ip00),dvzz(ip00),t(ip00)
      dimension sn(ip31*ndir)
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
!
      if(komsst.eq.0) then
!-----------------------------------------------------
!       modele de base Wilcox ou de Menter
!----------------------------------------------------
!
       do k=k1,k2m1
        do j=j1,j2m1
         n=ind(i1-1,j,k)
         do i=i1,i2m1
           n=n+nci
           mut(n)=v(n,1)*v(n,6)/v(n,7)
          enddo
         enddo
        enddo
!
      elseif(komsst.eq.1) then
!---------------------------------------------------
!       Modele Menter SST
!---------------------------------------------------
        betae2=betae/2.
        a1    =sqrt(betae)
!       Calcul du gradient de la vitesse. Les tableaux ont ete reutilises
!       dans la phase implicite.
!
        call teq_gradv( &
             l, &
             sn, &
             vol,v, &
             t , &
             dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
             cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)   !cEG-TOM
!
        do k=k1,k2m1
         do j=j1,j2m1
          n=ind(i1-1,j,k)
          do i=i1,i2m1
           n=n+nci
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
           mut(n)=v(n,6)/max(omeg,rota*f2/a1)
          enddo
         enddo
        enddo
      endif
!
      return
      end subroutine
end module
