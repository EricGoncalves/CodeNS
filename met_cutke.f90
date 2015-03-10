module mod_met_cutke
implicit none
contains
      subroutine met_cutke(l,s)
!
!***********************************************************************
!
!     Limiteurs inferieurs sur k et epsilon.
!     k devient "keps" et epsilon devient "epse" si k ou epsilon
!     est inferieur a la borne
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
implicit none
integer :: inc
integer :: indc
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: s
integer :: i1
integer :: i1m1
integer :: i2
integer :: i2m1
integer :: j1
integer :: j1m1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k1m1
integer :: k2
integer :: k2m1
integer :: n
integer :: n0c
integer :: nci
integer :: nid
integer :: nijd
integer :: njd
!
!-----------------------------------------------------------------------
!
      dimension s(ip11,ip60)
!
      indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
!
      n0c=npc(l)
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
      nci = inc(1,0,0)
!
      do k=k1,k2m1
       do j=j1,j2m1
        n=indc(i1m1,j,k)
        do i=i1,i2
         n=n+nci
         if(s(n,6).lt.epsk) then
          s(n,6)=epsk
          s(n,7)=epse
         endif
         if(s(n,7).lt.epse) then
          s(n,6)=epsk
          s(n,7)=epse
         endif
        enddo
       enddo
      enddo
!
      return
      end
end module
