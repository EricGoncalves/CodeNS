module mod_met_rfvc
  implicit none
contains
  subroutine met_rfvc(t,ncbd,mnc)
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use mod_mpi
    implicit none
    integer          ::          m,        mb,        mc,        mf,       mfb
    integer          ::  mnc(ip43),        mt,        nc,ncbd(ip41),        nd
    double precision :: t(ip11,ip60)
!    double precision,allocatable :: buff(:,:,:)
!
!-----------------------------------------------------------------------
!
!
!mt=0
!    do mf=1,nbd
!       mfb=lbd(mf)
!       mt=max(mt,mmb(mfb))
!    enddo
!allocate(buff(2,mt,nbd)

    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
!
      print*,"met_rfvc",rank,mfb,mt
       
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
!!!          buff(:,m,mf)=t(nc,6:7)
!!!!
!!!       enddo
!!!      MPI_trans(buff(:,1:mt,mf),buff(:,1:mt,mf),mf,mf)
!!!!
!!!    enddo

!!!    do mf=1,nbd
!!!!
!!!       mfb=lbd(mf)
!!!       mt=mmb(mfb)
!!!!
!!!       do m=1,mt
!!!          mc=mpc(mfb)+m
!!!          nc=mnc(mc)

!!!          mb=mpb(mfb)+m
!!!          nd=ncbd(mb)
!!!!
!!!!     definition des variables aux points fictifs
!!!!
!!!!          t(nd,6)=t(nc,6)
!!!!          t(nd,7)=t(nc,7)
!!!          t(nd,6:7)=buff(:,m,mf)
!
       enddo
!
    enddo
call barrier
stop
!
    return
  end subroutine met_rfvc
end module mod_met_rfvc
