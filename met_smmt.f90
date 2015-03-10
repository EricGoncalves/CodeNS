module mod_met_smmt
implicit none
contains
      subroutine met_smmt( &
                 l,ncyc, &
                 v,mu,mut,dist,mnpar,ncin, &
                 txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                 tprod,cfke, &
                 gkgo,frac, &
                 qcxts5,qcyts6)
!
!***********************************************************************
!
!_H   DATE_C :  Eric Goncalves / SINUMEF
!_H
!     ACT
!_A   Calcul du second membre des equations pour k-omega
!_A   Modele de Wilcox et Menter. Dans ce cas, la fonction de raccord
!_A   des modeles est aussi calculee.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_L    tprod      : arg real(ip00     )  ; production de k
!
!_O    frac       : arg real(ip12     )  ; fonction de raccord des modeles
!
!_/    txxf5x     : arg real(ip12     )  ; comp x grad(k) puis
!_/                                        tenseur visqueux
!_/    txyf5y     : arg real(ip12     )  ; comp y grad(k) puis
!_/                                        tenseur visqueux
!_/    txzf5z     : arg real(ip12     )  ; comp z grad(k) puis
!_/                                        tenseur visqueux
!_/    tyyf6x     : arg real(ip12     )  ; comp x grad(e) puis
!_/                                        tenseur visqueux
!_/    tyzf6y     : arg real(ip12     )  ; comp y grad(e) puis
!_/                                        tenseur visqueux
!_/    tzzf6z     : arg real(ip12     )  ; comp z grad(e) puis
!_/                                        tenseur visqueux
!_/    qcxts5     : arg real(ip12     )  ; terme source equation pour k puis
!_/                                        vecteur flux de chaleur
!_/    qcyts6     : arg real(ip12     )  ; terme source seconde equation puis
!
!_L    gkgo       : arg real(ip00     )  ; grad(k) * grad (omega)
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
use mod_met_mtcorf1
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
integer :: ncyc
double precision :: v
double precision :: dist
integer :: mnpar
integer :: ncin
double precision :: txxf5x
double precision :: txyf5y
double precision :: txzf5z
double precision :: tyyf6x
double precision :: tyzf6y
double precision :: tzzf6z
double precision :: tprod
double precision :: cfke
double precision :: gkgo
double precision :: frac
double precision :: qcxts5
double precision :: qcyts6
double precision :: coef1
double precision :: coef2
double precision :: coef3
double precision :: coef4
double precision :: dist2
double precision :: dkomg
double precision :: exp2x
double precision :: fm1
double precision :: fm2
integer :: i1
integer :: i1m1
integer :: i2
integer :: i2m1
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
integer :: n0c
integer :: nci
integer :: ncj
integer :: nck
integer :: nid
integer :: nijd
integer :: njd
double precision :: omeg
double precision :: sl3
double precision :: w4sig2
double precision :: wgam
double precision :: wsig
double precision :: xgam1
double precision :: xgam2
double precision :: zeta
!
!-----------------------------------------------------------------------
!
      real mu,mut
      dimension v(ip11,ip60)
      dimension mut(ip12),mu(ip12),dist(ip12),mnpar(ip12), &
                txxf5x(ip12),txyf5y (ip12),txzf5z(ip12), &
                tyyf6x(ip12),tyzf6y (ip12),tzzf6z(ip12), &
                qcxts5(ip12),qcyts6 (ip12),frac(ip12)
      dimension tprod(ip00),gkgo(ip00)
      dimension ncin(ip41),cfke(ip13)
!
      indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
!
      xgam1=beta1/betae-sigma1*okappa**2/sqrt(betae)
      xgam2=beta2/betae-sigma2*okappa**2/sqrt(betae)
      w4sig2=4.*wsig2
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
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
!
      nci  = inc(1,0,0)
      ncj  = inc(0,1,0)
      nck  = inc(0,0,1)
!     ----------------------------------------------------------
!com  calcul de la fonction de raccord entre modeles
!
!                       n : pointeur cellule tableaux tous domaines
!                       m : pointeur cellule tableaux un   domaine
!
      if(kfracom.eq.0) then
!
!       modele original de Wilcox
!
        do k=k1,k2m1
         do j=j1,j2m1
          n=indc(i1m1,j,k)
          do i=i1,i2m1
           n=n+nci
           m=n-n0c
           frac(n)=1.
           gkgo(m)=0.
          enddo
         enddo
        enddo
      else
!
!        modele de Menter
!
        do k=k1,k2m1
         do j=j1,j2m1
          n=indc(i1m1,j,k)
          do i=i1,i2m1
           n=n+nci
           m=n-n0c
           gkgo(m)=tyyf6x(n)*txxf5x(n)+tyzf6y(n)*txyf5y(n)+ &
                   tzzf6z(n)*txzf5z(n)
           omeg=v(n,7)/v(n,1)
           coef1=sqrt(v(n,6)/v(n,1))/(betae*omeg*dist(n))
           dist2=dist(n)**2
           coef2=500.*mu(n)/(v(n,7)*dist2)
           coef4=wsig2*v(n,1)*gkgo(m)/omeg
           dkomg=max(coef4,1.e-20)
           coef3=w4sig2*v(n,6)/(dist2*dkomg)
           zeta =min( max(coef1,coef2),coef3)
           if(zeta.le.2.5) then
            exp2x=exp(2.*zeta**4)
            frac(n)=(exp2x-1.)/(exp2x+1.)
           else
            frac(n)=1.
           endif
          enddo
         enddo
        enddo
!
!        verification que la fonction de raccord reste a 1 dans la region
!        de paroi
!
         call met_mtcorf1( &
                 l,ncin, &
                 dist,mnpar,frac)
!
      endif
!     --------------------------------------------------------------
!com  calcul des termes sources et du rayon spectral pour implicitation
!
      do k=k1,k2m1
       do j=j1,j2m1
         n=indc(i1m1,j,k)
         do i=i1,i2m1
          n=n+nci
          m=n-n0c
          fm1=frac(n)
          fm2=1.-fm1
          wgam=fm1*xgam1+fm2*xgam2
          wsig=fm1*wsig1+fm2*wsig2
          beta=fm1*beta1+fm2*beta2
          omeg=v(n,7)/v(n,1)
          sl3=2.*wsig*v(n,1)*gkgo(m)/omeg
          qcxts5(n)=tprod(m) -betae*v(n,6)*omeg
          qcyts6(n)=tprod(m)*wgam*v(n,1)/mut(n) &
                  - beta*v(n,7)*omeg + sl3
          cfke(n)=0.18*omeg
         enddo
        enddo
      enddo
!
      return
      end
end module
