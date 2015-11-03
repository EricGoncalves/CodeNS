module mod_met_fludko
  implicit none
contains
  subroutine met_fludko( &
       l,mu,mut, &
       fd5x,fd5y,fd5z,fd6x,fd6y,fd6z)
!
!***********************************************************************
!
!_DA  DATE_C :  avril 2002    - Eric Goncalves / SINUMEF
!
!     ACT
!_A   Calcul des densites de flux dissipatifs a partir de grad(k) et grad(w)
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    sigmk      : com real             ; constante sigma_k
!_I    sigmw      : com real             ; constante sigma_omega
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
    implicit none
    integer          ::    i,  i1,  i2,i2m1,ind1
    integer          :: ind2,   j,  j1,  j2,j2m1
    integer          ::    k,  k1,  k2,k2m1,   l
    integer          ::    n,  n0, nid,nijd, njd
    double precision :: fd5x,fd5y,fd5z,fd6x,fd6y
    double precision :: fd6z,  mu, mut
!
!-----------------------------------------------------------------------
!
    dimension mut(ip12),mu(ip12), &
         fd5x(ip12),fd5y(ip12),fd5z(ip12), &
         fd6x(ip12),fd6y(ip12),fd6z(ip12)
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
    ind1=indc(i1  ,j1  ,k1  )
    ind2=indc(i2m1,j2m1,k2m1)
    do n=ind1,ind2
       fd5x(n)=(mu(n)+mut(n)*sigmk)*fd5x(n)
       fd5y(n)=(mu(n)+mut(n)*sigmk)*fd5y(n)
       fd5z(n)=(mu(n)+mut(n)*sigmk)*fd5z(n)
!
       fd6x(n)=(mu(n)+mut(n)*sigmw)*fd6x(n)
       fd6y(n)=(mu(n)+mut(n)*sigmw)*fd6y(n)
       fd6z(n)=(mu(n)+mut(n)*sigmw)*fd6z(n)
    enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
  end subroutine met_fludko
end module mod_met_fludko
