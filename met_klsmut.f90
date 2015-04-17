module mod_met_klsmut
  implicit none
contains
  subroutine met_klsmut( &
       l, &
       sn,vol,t, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       dist,v,mu,mut, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DA  DATE_C :  juin 2000 - Eric Goncalves
!
!     ACT
!_A   Calcul de la viscosite turbulente a partir de k et l
!_A   Modele de Smith - modele avec correction SST
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use mod_teq_gradv
    implicit none
    integer          ::    i,  i1,  i2,i2m1,  id
    integer          ::    j,  j1,  j2,j2m1,  jd
    integer          ::    k,  k1,  k2,k2m1,  kd
    integer          ::    l,   m,   n,  n0, nci
    integer          ::  ncj, nck, nid,nijd, njd
    double precision ::            a1,           c1,          c14,           c2,          c22
    double precision ::   cmui1(ip21),  cmui2(ip21),  cmuj1(ip21),  cmuj2(ip21),  cmuk1(ip21)
    double precision ::   cmuk2(ip21),        coef1,        coef2,   dist(ip12),   dvxx(ip00)
    double precision ::    dvxy(ip00),   dvxz(ip00),   dvyx(ip00),   dvyy(ip00),   dvyz(ip00)
    double precision ::    dvzx(ip00),   dvzy(ip00),   dvzz(ip00),        exp2x,           f1
    double precision ::            f2,          fmu,     mu(ip12),    mut(ip12),        r2sb1
    double precision ::          rack,         rota,sn(ip31*ndir),          ss2,      t(ip00)
    double precision ::  v(ip11,ip60),    vol(ip11),          xl2,          xxi,         xxi2
    double precision ::          zeta
!
!-----------------------------------------------------------------------
!
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
!
    r2sb1 =sqrt(2.)/cklb1**(1./3.)
    c14=c1**4
    c22=c2**2
!
!     correction SST
    a1=sqrt(cmukl) !a1=0.3
!
!       Calcul du gradient de la vitesse. Les tableaux ont ete reutilises
!       dans la phase implicite.
!
    call teq_gradv( &
         l, &
         sn, &
         vol,v, &
         t , &
         dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          n=ind(i1-1,j,k)
!$OMP SIMD
          do i=i1,i2m1
             n=n+nci
             m=n-n0
             rack=sqrt(v(n,6)/v(n,1))
             xxi=r2sb1*rack*v(n,7)/mu(n)
             xxi2=xxi*xxi
             xl2=(v(n,7)/(xkappa*v(n,1)*dist(n)))**2
             f1=exp(-50.*xl2)
             ss2=(c22+xxi2)*xxi2
             fmu=((c14*f1+ss2)/(c14+ss2))**0.25
!
             coef1=2.*r2sb1*v(n,7)/(cmukl*dist(n)*v(n,1))
             coef2=500.*r2sb1*mu(n)*v(n,7)/(rack*v(n,1)**2*dist(n)**2)
             rota = sqrt( (dvzy(m)-dvyz(m))**2 &
                  +(dvxz(m)-dvzx(m))**2 &
                  +(dvyx(m)-dvxy(m))**2)
             zeta=max(coef1,coef2)
             exp2x=exp(min(2.*zeta**2,25.))
             f2=(exp2x-1.)/(exp2x+1.)
             a1=0.3
             mut(n)=v(n,6)/max(rack*v(n,1)/(fmu*v(n,7)*r2sb1),rota*f2/a1)
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
  end subroutine met_klsmut
end module mod_met_klsmut
