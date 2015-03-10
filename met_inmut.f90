      subroutine met_inmut(l,mu,mut)
!
!***********************************************************************
!
!        initalisation mut = mu/10
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
implicit none
integer :: indc
integer :: i
integer :: j
integer :: k
integer :: l
integer :: i1
integer :: i2
integer :: i2m1
integer :: ind1
integer :: ind2
integer :: j1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k2
integer :: k2m1
integer :: n
integer :: n0c
integer :: nid
integer :: nijd
integer :: njd
!
!-----------------------------------------------------------------------
!
      real mu,mut
      dimension mu(ip12),mut(ip12)
!
      indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
!
      n0c=npc(l)
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
      do k=k1,k2m1
       do j=j1,j2m1
        ind1=indc(i1,j,k)
        ind2=indc(i2m1,j,k)
        do n=ind1,ind2
         mut(n)=mu(n)*1.e-1
        enddo
       enddo
      enddo
!
      return
      end
