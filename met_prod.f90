module mod_met_prod
  implicit none
contains
  subroutine met_prod( &
       l,ncyc, &
       mut,s, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       tprod)
!
!***********************************************************************
!
!     ACT  calcul de la production de k
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::    i,  i1,i1m1,  i2,i2m1
    integer          :: ind1,ind2,   j,  j1,j1m1
    integer          ::   j2,j2m1,   k,  k1,k1m1
    integer          ::   k2,k2m1,   l,   m,   n
    integer          ::  n0c,ncyc, nid,nijd, njd
    double precision ::   ds3, dvxx, dvxy, dvxz, dvyx
    double precision ::  dvyy, dvyz, dvzx, dvzy, dvzz
    double precision ::   mut,    s,torxx,torxy,torxz
    double precision :: toryy,toryz,torzz,tprod
!
!-----------------------------------------------------------------------
!
    dimension s(ip11,ip60)
    dimension mut(ip12)
    dimension tprod(ip00), &
         dvxx(ip00),dvxy(ip00),dvxz(ip00), &
         dvyx(ip00),dvyy(ip00),dvyz(ip00), &
         dvzx(ip00),dvzy(ip00),dvzz(ip00)
!

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
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    ds3=2./3.
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             torxx=ds3*(-s(n,6)-mut(n)*(-2*dvxx(m)+dvyy(m)+dvzz(m)))
             toryy=ds3*(-s(n,6)-mut(n)*(-2*dvyy(m)+dvxx(m)+dvzz(m)))
             torzz=ds3*(-s(n,6)-mut(n)*(-2*dvzz(m)+dvyy(m)+dvxx(m)))
             torxy=             mut(n)*(dvxy(m)+dvyx(m))
             torxz=             mut(n)*(dvxz(m)+dvzx(m))
             toryz=             mut(n)*(dvzy(m)+dvyz(m))
             tprod(m)=dvxx(m)*torxx +dvyy(m)*toryy +dvzz(m)*torzz &
                  +torxy*(dvxy(m)+dvyx(m)) &
                  +torxz*(dvxz(m)+dvzx(m)) &
                  +toryz*(dvzy(m)+dvyz(m))
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
  end subroutine met_prod
end module mod_met_prod
