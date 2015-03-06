      subroutine met_cut(l,s)
!
!***********************************************************************
!
!     ACT
!      Application des limiteurs sur les variables turbulentes.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
!
!-----------------------------------------------------------------------
!
      real test6,test7,ts
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
        n=indc(i1-1,j,k)
        do i=i1,i2m1
         n=n+nci
         test6=sign(1.,s(n,6)-epsk)
         test7=sign(1.,s(n,7)-epse)
         ts=sign(1.,test6+test7-1.5)
         s(n,6)=0.5*((1.+ts)*s(n,6)+(1.-ts)*epsk)
         s(n,7)=0.5*((1.+ts)*s(n,7)+(1.-ts)*epse)
        enddo
       enddo
      enddo
!
      return
      end
