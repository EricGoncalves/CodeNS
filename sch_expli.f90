module mod_sch_expli
  implicit none
contains
  subroutine sch_expli(lm,u,v,dt,vol)
!
!***********************************************************************
!
!_DA  DATE_C : avril 2003 - Eric Goncalves / LEGI
!
!     ACT
!       Avance en temps explicite
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::   i1,i1m1,i1p1,  i2
    integer          :: i2m1,ind1,ind2,   j,  j1
    integer          :: j1m1,j1p1,  j2,j2m1,   k
    integer          ::   k1,k1m1,k1p1,  k2,k2m1
    integer          ::   lm,   n, n0c, nid,nijd
    integer          ::  njd
    double precision ::     dt(ip11),         dtv,u(ip11,ip60),v(ip11,ip60),   vol(ip11)
!
!-----------------------------------------------------------------------
!
!
    n0c=npc(lm)
    i1=ii1(lm)
    i2=ii2(lm)
    j1=jj1(lm)
    j2=jj2(lm)
    k1=kk1(lm)
    k2=kk2(lm)
!
    nid = id2(lm)-id1(lm)+1
    njd = jd2(lm)-jd1(lm)+1
    nijd = nid*njd
!
    i1p1=i1+1
    j1p1=j1+1
    k1p1=k1+1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             dtv=dt(n)/vol(n)
             v(n,1)=v(n,1)+dtv*u(n,1)
             v(n,2)=v(n,2)+dtv*u(n,2)
             v(n,3)=v(n,3)+dtv*u(n,3)
             v(n,4)=v(n,4)+dtv*u(n,4)
             v(n,5)=v(n,5)+dtv*u(n,5)
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
  end subroutine sch_expli
end module mod_sch_expli
