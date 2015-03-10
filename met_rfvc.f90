      subroutine met_rfvc(t,ncbd,mnc)
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
   use boundary
!
!-----------------------------------------------------------------------
!
      dimension t(ip11,ip60)
      dimension mnc(ip43),ncbd(ip41)
!
      do mf=1,nbd
!
         mfb=lbd(mf)
         mt=mmb(mfb)
!
!DEC$ IVDEP
         do m=1,mt
            mc=mpc(mfb)+m
            nc=mnc(mc)
            mb=mpb(mfb)+m
            nd=ncbd(mb)
!
!     definition des variables aux points fictifs
!
            t(nd,6)=t(nc,6)
            t(nd,7)=t(nc,7)
!
         enddo
!
      enddo
!
      return
      end
