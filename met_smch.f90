module mod_met_smch
implicit none
contains
      subroutine met_smch( &
                 l, &
                 v,mu,dist,mnpar,utau, &
                 tprod, &
                 qcxts5,qcyts6)
!
!***********************************************************************
!
!     ACT
!_A   Modele de Chien
!_A   Calcul du second membre des equations pour k-epsilon
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_/    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_L    tprod      : arg real(ip00     )  ; production de k
!_L    dvxx       : arg real(ip00     )  ; grad(V)  vx,x
!_L    dvxy       : arg real(ip00     )  ; grad(V)  vx,y
!_L    dvxz       : arg real(ip00     )  ; grad(V)  vx,z
!_L    dvyx       : arg real(ip00     )  ; grad(V)  vy,x
!_L    dvyy       : arg real(ip00     )  ; grad(V)  vy,y
!_L    dvyz       : arg real(ip00     )  ; grad(V)  vy,z
!_L    dvzx       : arg real(ip00     )  ; grad(V)  vz,x
!_L    dvzy       : arg real(ip00     )  ; grad(V)  vz,y
!_L    dvzz       : arg real(ip00     )  ; grad(V)  vz,z
!
!     OUT
!_O
!     I/O
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
!     LOC
!_L    t          : arg real(ip00     )  ; variable de travail
!_L    dtdx       : arg real(ip00     )  ; grad(t)  t,x
!_L    dtdy       : arg real(ip00     )  ; grad(t)  t,y
!_L    dtdz       : arg real(ip00     )  ; grad(t)  t,z
!_L    bark       : arg real(ip00     )  ; terme faible reynolds equation k
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
      use chainecarac
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
double precision :: v
double precision :: dist
integer :: mnpar
double precision :: utau
double precision :: tprod
double precision :: qcxts5
double precision :: qcyts6
double precision :: dist2
double precision :: epssk
double precision :: f2
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
integer :: mp
integer :: n
integer :: n0c
integer :: nci
integer :: ncj
integer :: nck
integer :: nid
integer :: nijd
integer :: njd
double precision :: rodst
double precision :: rtur
double precision :: sch1
double precision :: sch2
double precision :: sch3
double precision :: xk
double precision :: yplus
!
!-----------------------------------------------------------------------
!
      real mu
      dimension v(ip11,ip60)
      dimension mu(ip12),dist(ip12),mnpar(ip12),utau(ip42), &
                qcxts5(ip12),qcyts6 (ip12)
      dimension tprod(ip00)
!
      indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
!
!     ----------------------------------------------------------
!
!     notation pour les constantes
!
!     c_epsilon_1 <-> cke1
!     c_epsilon_2 <-> cke2
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
      ncj  = inc(0,1,0)
      nck  = inc(0,0,1)
      do k=k1,k2m1
         do j=j1,j2m1
            n=indc(i1m1,j,k)
            do i=i1,i2m1
               n    =n+nci
               m    =n-n0c
               mp   =mnpar(n)
               if(mp.ge.1) then
                 yplus=max(abs(utau(mp)),utaumin)*dist(n)*v(n,1)/mu(n)
               else
                 yplus=1000.
               end if
               xk   =v(n,6)/v(n,1)
               dist2=dist(n)**2
               rodst=v(n,1)*dist2
               qcxts5(n)=tprod(m) - v(n,7) - 2.*mu(n)*v(n,6)/rodst
!
               epssk=v(n,7)/v(n,6)
               rtur =v(n,6)/(mu(n)*epssk)
               f2   =1.-.22222222*exp(-rtur*rtur/36.)
               sch1 =cke1*epssk*tprod(m)
               sch2 =cke2*f2*epssk*v(n,7)
               sch3 =2.*mu(n)*(v(n,7)/rodst)*exp(-yplus/2.)
               qcyts6(n)=sch1 - sch2 - sch3
            enddo
         enddo
      enddo
!
      return
      end subroutine
end module
