module mod_sch_dual
  implicit none
contains
  subroutine sch_dual( &
       lm,u,v,icycle, &
       vol,ptdual)
!
!***********************************************************************
!
!_DA  DATE_C : octobre 2001 : Eric Goncalves - SINUMEF
!
!     ACT
!_A    Calcul du residu instationnaire R* pour ordre 2
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
    integer          ::     lm,     n,   n0c,   nid,  nijd
    integer          ::    njd
    double precision ::                c0,              dti,             fact,ptdual(ip11,ip60),     u(ip11,ip60)
    double precision ::      v(ip11,ip60),        vol(ip11)
!
!-----------------------------------------------------------------------
!
!

!
    fact=1.5
    if(icycle.le.1) fact=1.
    dti=1./dt1min
!
    n0c=npc(lm)
    i1=ii1(lm)
    i2=ii2(lm)
    j1=jj1(lm)
    j2=jj2(lm)
    k1=kk1(lm)
    k2=kk2(lm)
    nid = id2(lm)-id1(lm)+1
    njd = jd2(lm)-jd1(lm)+1
    nijd = nid*njd
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1  ,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             c0=vol(n)*dti
             u(n,1)=u(n,1)+c0*(fact*v(n,1)+ ptdual(n,1))
             u(n,2)=u(n,2)+c0*(fact*v(n,2)+ ptdual(n,2))
             u(n,3)=u(n,3)+c0*(fact*v(n,3)+ ptdual(n,3))
             u(n,4)=u(n,4)+c0*(fact*v(n,4)+ ptdual(n,4))
             u(n,5)=u(n,5)+c0*(fact*v(n,5)+ ptdual(n,5))
          enddo
       enddo
    enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
    end function indc
  end subroutine sch_dual
end module mod_sch_dual
