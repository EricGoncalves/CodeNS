module mod_met_klmut
  implicit none
contains
  subroutine met_klmut( &
       l, &
       v,mu,mut,dist)
!
!***********************************************************************
!
!     ACT
!_A   Calcul de la viscosite turbulente a partir de k et l
!_A   Modele de Smith
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use proprieteflu
    implicit none
    integer          ::    i,  i1,  i2,i2m1, ncj
    integer          ::    j,  j1,  j2,j2m1
    integer          ::    k,  k1,  k2,k2m1
    integer          ::    l,   n,  n0, nci
    integer          ::  nck, nid,nijd, njd
    double precision ::           c1,         c14,          c2,         c22,  dist(ip12)
    double precision ::           f1,         fmu,    mu(ip12),   mut(ip12),       r2sb1
    double precision ::         rack,         ss2,v(ip11,ip60),         xl2,         xxi
    double precision ::         xxi2
!
!-----------------------------------------------------------------------
!
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
!     notation pour les constantes
!
!     B1 <-> cklb1
!     E2 <-> ckle2
!     sigma_k <-> sigmak
!     sigma_l <-> sigmal
!
    c1=25.5
    c2=2.
    r2sb1 =sqrt(2.)/cklb1**(1./3.)
    c14=c1**4
    c22=c2**2
!
    do k=k1,k2m1
       do j=j1,j2m1
          n=ind(i1-1,j,k)
          do i=i1,i2m1
             n=n+nci
             rack=sqrt(v(n,6)/v(n,1))
             xxi=r2sb1*rack*v(n,7)/mu(n)
             xxi2=xxi*xxi
             xl2=(v(n,7)/(xkappa*v(n,1)*dist(n)))**2
             f1=exp(-50.*xl2)
             ss2=(c22+xxi2)*xxi2
             fmu=((c14*f1+ss2)/(c14+ss2))**0.25
             mut(n)=fmu*mu(n)*xxi
          enddo
       enddo
    enddo
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine met_klmut
end module mod_met_klmut
