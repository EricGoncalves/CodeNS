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
    integer          ::   i1,  i2,i2m1, img
    integer          :: ind1,ind2,   j,  j1,  j2
    integer          :: j2m1,   k,  k1,  k2,k2m1
    integer          ::    l,  lm,   m, n0c,  nc
    integer          ::  nid,nijd, njd
    double precision :: dc(ip11,ip60), v(ip11,ip60)
!
!-----------------------------------------------------------------------
!
!

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
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
    end function indc
  end subroutine smg_upc
end module mod_smg_upc
