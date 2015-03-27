module mod_met_dual
  implicit none
contains
  subroutine met_dual( &
       icycle,l,u,v, &
       vol,ptdual)
!
!***********************************************************************
!
!_DA  DATE_C : octobre 2001 : Eric GONCALVES - SINUMEF
!
!     ACT
!_A    Calcul du residu instationnaire R* pour les equations
!-A    de transport de la turbulence - ordre 2
!
!-----------------------------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use schemanum
    implicit none
    integer          ::      i,    i1,    i2,  i2m1,icycle
    integer          ::   ind1,  ind2,     j,    j1,    j2
    integer          ::   j2m1,     k,    k1,    k2,  k2m1
    integer          ::      l,     n,   n0c,   nid,  nijd
    integer          ::    njd
    double precision ::                c0,              dti,             fact,ptdual(ip11,ip60),     u(ip11,ip60)
    double precision ::      v(ip11,ip60),        vol(ip11)
!
!-----------------------------------------------------------------------
!
!

!
    fact = 1.5
    if(icycle.eq.1) fact = 1.
    dti = 1./dt1min
!
    n0c  = npc(l)
    i1   = ii1(l)
    i2   = ii2(l)
    j1   = jj1(l)
    j2   = jj2(l)
    k1   = kk1(l)
    k2   = kk2(l)
    nid  = id2(l)-id1(l)+1
    njd  = jd2(l)-jd1(l)+1
    nijd = nid*njd
    i2m1 = i2-1
    j2m1 = j2-1
    k2m1 = k2-1
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
          do n=ind1,ind2
             c0=vol(n)*dti
             u(n,6)=u(n,6) + c0*(fact*v(n,6)+ptdual(n,6))
             u(n,7)=u(n,7) + c0*(fact*v(n,7)+ptdual(n,7))
          enddo
       enddo
    enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
  end subroutine met_dual
end module mod_met_dual
