module mod_met_fludmt
implicit none
contains
      subroutine met_fludmt( &
                 l, &
                 s,mu,mut, &
                 frac, &
                 fd5x,fd5y,fd5z,fd6x,fd6y,fd6z)
!
!***********************************************************************
!
!     ACT
!_A   Calcul des densites de flux dissipatifs a partir de grad(k) et grad(e)
!_A   pour "sigma_k" et "sigma_e" constants.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    s          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    qczfrc     : arg real(ip12      ) ; fonction de raccord de Menter (F1)
!_I    cmu        : com real             ; constante C_mu
!_I    cke1       : com real             ; constante C_epsilon-1
!_I    cke2       : com real             ; constante C_epsilon-2
!_I    alfak      : com real             ; constante sigma_k
!_I    alfak      : com real             ; constante sigma_epsilon
!
!     I/O
!_/    fd5x       : arg real(ip12     )  ; comp x grad(k) puis flux diss.
!_/    fd5y       : arg real(ip12     )  ; comp y grad(k) puis flux diss.
!_/    fd5z       : arg real(ip12     )  ; comp z grad(k) puis flux diss.
!_/    fd6x       : arg real(ip12     )  ; comp x grad(e) puis flux diss.
!_/    fd6y       : arg real(ip12     )  ; comp y grad(e) puis flux diss.
!_/    fd6z       : arg real(ip12     )  ; comp z grad(e) puis flux diss.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
      use chainecarac
implicit none
integer :: indc
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: s
double precision :: frac
double precision :: fd5x
double precision :: fd5y
double precision :: fd5z
double precision :: fd6x
double precision :: fd6y
double precision :: fd6z
double precision :: fm1
double precision :: fm2
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
double precision :: sigme
!
!-----------------------------------------------------------------------
!
      real mu,mut,fd5x0
      dimension s(ip11,ip60)
      dimension mut(ip12),mu(ip12),frac(ip12), &
                fd5x(ip12),fd5y(ip12),fd5z(ip12), &
                fd6x(ip12),fd6y(ip12),fd6z(ip12)
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
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
!
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd= nid*njd
!
      ind1=indc(i1  ,j1  ,k1  )
      ind2=indc(i2m1,j2m1,k2m1)
!
      if(kfracom.eq.0) then
!
!       modele de Wilcox
!
        do n=ind1,ind2
           m=n-n0
           fd5x0=fd5x(n)
           fd5x(n)=(mu(n)+mut(n)*sigme1)*fd5x(n)
           fd5y(n)=(mu(n)+mut(n)*sigme1)*fd5y(n)
           fd5z(n)=(mu(n)+mut(n)*sigme1)*fd5z(n)
!
           fd6x(n)=(mu(n)+mut(n)*sigma1)*fd6x(n)
           fd6y(n)=(mu(n)+mut(n)*sigma1)*fd6y(n)
           fd6z(n)=(mu(n)+mut(n)*sigma1)*fd6z(n)
        enddo
      else
!
!       Modele de Menter
!
        do n=ind1,ind2
           m=n-n0
           fm1=frac(n)
           fm2=1.-fm1
           sigme=fm1*sigme1+fm2*sigme2
           sigma=fm1*sigma1+fm2*sigma2
           fd5x(n)=(mu(n)+mut(n)*sigme)*fd5x(n)
           fd5y(n)=(mu(n)+mut(n)*sigme)*fd5y(n)
           fd5z(n)=(mu(n)+mut(n)*sigme)*fd5z(n)
!
           fd6x(n)=(mu(n)+mut(n)*sigma)*fd6x(n)
           fd6y(n)=(mu(n)+mut(n)*sigma)*fd6y(n)
           fd6z(n)=(mu(n)+mut(n)*sigma)*fd6z(n)
        enddo
      endif
!
      return
      end subroutine
end module
