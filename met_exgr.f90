module mod_met_exgr
  implicit none
contains
  subroutine met_exgr( &
       l, &
       ds1,ds2,ds3)
!***********************************************************************
!
!     ACT
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::    i,  i1,i1m1,i1p1,  i2
    integer          :: i2m1,i2p1,  id,iinc,  is
    integer          ::    j,  j1,j1m1,j1p1,  j2
    integer          :: j2m1,j2p1,  jd,jinc,  js
    integer          ::    k,  k1,k1m1,k1p1,  k2
    integer          :: k2m1,k2p1,  kd,kinc,  ks
    integer          ::    l,   m,   n, n0c, nci
    integer          ::  ncj, nck, nid,nijd, njd
    double precision :: ds1(ip00),ds2(ip00),ds3(ip00),      ex1,      ex2
!
!-----------------------------------------------------------------------
!
!
!


    ex1=1.
    ex2=0.
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
!
    i2p1=i2+1
    j2p1=j2+1
    k2p1=k2+1
!
    i1p1=i1+1
    j1p1=j1+1
    k1p1=k1+1
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    nci  = inc(1,0,0)
    ncj  = inc(0,1,0)
    nck  = inc(0,0,1)
!
! extrapolation des variables
!
    iinc=i2-i1m1
    jinc=j2-j1m1
    kinc=k2-k1m1
!
    ks=-1
    do k=k1m1,k2,kinc
       ks=-ks
       do j=j1,j2m1
!$OMP SIMD
          do i=i1,i2m1
             n=indc(i,j,k)
             m=n-n0c
             ds1(m)=ex1*ds1(m+ks*nck)+ex2*ds1(m+2*ks*nck)
             ds2(m)=ex1*ds2(m+ks*nck)+ex2*ds2(m+2*ks*nck)
             ds3(m)=ex1*ds3(m+ks*nck)+ex2*ds3(m+2*ks*nck)
          enddo
       enddo
    enddo
!
    do k=k1m1,k2
       do j=j1,j2m1
          is=-1
          do i=i1m1,i2,iinc
             is=-is
             n=indc(i,j,k)
             m=n-n0c
             ds1(m)=ex1*ds1(m+is*nci)+ex2*ds1(m+2*is*nci)
             ds2(m)=ex1*ds2(m+is*nci)+ex2*ds2(m+2*is*nci)
             ds3(m)=ex1*ds3(m+is*nci)+ex2*ds3(m+2*is*nci)
          enddo
       enddo
    enddo
!
    do k=k1m1,k2
       js=-1
       do j=j1m1,j2,jinc
          js=-js
!$OMP SIMD
          do i=i1m1,i2
             n=indc(i,j,k)
             m=n-n0c
             ds1(m)=ex1*ds1(m+js*ncj)+ex2*ds1(m+2*js*ncj)
             ds2(m)=ex1*ds2(m+js*ncj)+ex2*ds2(m+2*js*ncj)
             ds3(m)=ex1*ds3(m+js*ncj)+ex2*ds3(m+2*js*ncj)
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
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine met_exgr
end module mod_met_exgr
