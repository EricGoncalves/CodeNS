module mod_met_brad
  implicit none
contains
  subroutine met_brad( &
       l, &
       mut,v, &
       dvxy,dvxz,dvyx,dvyz,dvzx,dvzy)
!
!***********************************************************************
!
!     ACT
!     Utilisation de la vorticite pour initaliser k et epsilon
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    implicit none
    integer          ::    i,  i1,i1m1,  i2,i2m1
    integer          :: ind1,ind2,   j,  j1,j1m1
    integer          ::   j2,j2m1,   k,  k1,k1m1
    integer          ::   k2,k2m1,   l,   m,   n
    integer          ::   n0, nid,nijd, njd
    double precision ::         cmu1,  dvxy(ip00),  dvxz(ip00),  dvyx(ip00),  dvyz(ip00)
    double precision ::   dvzx(ip00),  dvzy(ip00),   mut(ip12),        omeg,v(ip11,ip60)
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
    cmu1=1./sqrt(cmu)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0
             omeg = sqrt((dvzy(m)-dvyz(m))**2 &
                  +(dvxz(m)-dvzx(m))**2 &
                  +(dvyx(m)-dvxy(m))**2)
             v(n,6)=cmu1*mut(n)*omeg
             v(n,7)=mut(n)*omeg*omeg
          enddo
       enddo
    enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
  end subroutine met_brad
end module mod_met_brad
