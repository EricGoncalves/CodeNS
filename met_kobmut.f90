module mod_met_kobmut
implicit none
contains
      subroutine met_kobmut( &
                 ncyc, &
                 l, &
                 v,mu,mut)
!
!***********************************************************************
!
!_DA  DATE_C :  novembre 1998 - AUTEUR :  Eric Goncalves
!
!     ACT
!_A   Calcul de la viscosite turbulente a partir de k et omega
!_A   Modele k-omega bas Reynolds de Wilcox
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
implicit none
integer :: inc
integer :: ind
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: ncyc
integer :: l
double precision :: v
double precision :: allfae
integer :: i1
integer :: i2
integer :: i2m1
integer :: idebug
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
double precision :: omegha
double precision :: rapk
double precision :: reytur
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
      idebug=0
!c
!c      modele bas Reynolds de Wilcox ou de Menter
!c
        do k=k1,k2m1
         do j=j1,j2m1
          n=ind(i1-1,j,k)
          do i=i1,i2m1
           n=n+nci
         omegha=v(n,7)/v(n,1)
         reytur=v(n,6)/(omegha*mu(n))
           rapk=(reytur/rrk)
         allfae=(allfae0+rapk)/(1.+rapk)
           mut(n)=allfae*v(n,6)/omegha
          enddo
         enddo
        enddo
!
      return
      end subroutine
end module
