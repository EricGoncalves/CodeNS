module mod_met_samut
implicit none
contains
      subroutine met_samut(l,v,mu,mut)
!
!***********************************************************************
!
!     ACT
!_A   modele de Spalart et Allmaras : calcul du mu_t
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!
!     OUT
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
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
double precision :: v
double precision :: cv13
double precision :: fv1
integer :: i1
integer :: i2
integer :: i2m1
integer :: j1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k2
integer :: k2m1
integer :: n
integer :: n0
integer :: nci
integer :: ncj
integer :: nck
integer :: nid
integer :: nijd
integer :: njd
double precision :: xkhi
double precision :: xkhi3
!
!-----------------------------------------------------------------------
!
      double precision mu,mut
      dimension mu(ip12),mut(ip12)
      dimension v(ip11,ip60)
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
      cv13=cv1**3
      do k=k1,k2m1
       do j=j1,j2m1
        n=ind(i1-1,j,k)
        do i=i1,i2m1
         n=n+nci
         xkhi=v(n,6)/mu(n)
         xkhi3=xkhi**3
         fv1=xkhi3/(xkhi3+cv13)
         mut(n)=v(n,6)*fv1
        enddo
       enddo
      enddo
!
      return
      end subroutine
end module
