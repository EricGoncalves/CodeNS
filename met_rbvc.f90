      subroutine met_rbvc( &
                 t, &
                 ncbd,ncin,mnc)
!
!***********************************************************************
!
!     ACT
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
      dimension ncbd(ip41)
      dimension ncin(ip41)
      dimension mnc(ip43)
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
            ndm=ncin(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
            t(nd,6) = 0.5*( t(ndm,6)+t(nc,6) )
            t(nd,7) = 0.5*( t(ndm,7)+t(nc,7) )
!
         enddo
!
      enddo
!
      return
      end
