module mod_smg_upc
implicit none
contains
      subroutine smg_upc( &
                        img, &
                        dc,v)
!
!***********************************************************************
!
!     ACT
!_A    Update the Corrections from the Coarse to Fine level.
!_A    Corrections are updated on the next finer level "img".
!
!***********************************************************************
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
integer :: img
double precision :: dc
double precision :: v
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
integer :: l
integer :: lm
integer :: m
integer :: n0c
integer :: nc
integer :: nid
integer :: nijd
integer :: njd
!
!-----------------------------------------------------------------------
!
      dimension dc(ip11,ip60),v(ip11,ip60)
!
      indc(i,j,k)=1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
!
      do l = 1,lzx
      lm = l+(img-1)*lz
!
      n0c=npc(lm)
      i1 =ii1(lm)
      i2 =ii2(lm)
      j1 =jj1(lm)
      j2 =jj2(lm)
      k1 =kk1(lm)
      k2 =kk2(lm)
!
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
      nid  = id2(lm)-id1(lm)+1
      njd  = jd2(lm)-jd1(lm)+1
      nijd = nid*njd
!
      do k = k1,k2m1
       do j = j1,j2m1
        ind1 = indc(i1,j,k)
        ind2 = indc(i2m1,j,k)
        do m = ind1,ind2
         nc= m+n0c
         v(nc,1) = v(nc,1) + dc(nc,1)
         v(nc,2) = v(nc,2) + dc(nc,2)
         v(nc,3) = v(nc,3) + dc(nc,3)
         v(nc,4) = v(nc,4) + dc(nc,4)
         v(nc,5) = v(nc,5) + dc(nc,5)
        enddo
       enddo
      enddo
!
      enddo
!
      return
      end subroutine
end module
