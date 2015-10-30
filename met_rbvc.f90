module mod_met_rbvc
  implicit none
contains
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
    use mod_mpi
    implicit none
    integer          ::          m,        mb,        mc,        mf,       mfb
    integer          ::  mnc(ip43),        mt,        nc,ncbd(ip41),ncin(ip41)
    integer          ::         nd,       ndm
    double precision :: t(ip11,ip60)
    double precision,allocatable :: buff(:,:,:,:)
    integer :: req(nbd,2),other,me,bcg_to_mf(num_bcl)
!
!-----------------------------------------------------------------------
!
!
    req=MPI_REQUEST_NULL
    mt=0
    do mf=1,nbd
       mfb=lbd(mf)
       mt=max(mt,mmb(mfb))
       me=bcl_to_bcg(mfb)
       bcg_to_mf(me)=mf
    enddo
    allocate(buff(2,mt,nbd,2))

    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
       other=ndcc(mfb)
       me=bcl_to_bcg(mfb)
!
!     we have to exchange the globally numbered me boundary with the owner of the globally numbered other boundary
!
       do m=1,mt
          mb=mpb(mfb)+m
          nc=ncin(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
          buff(1,m,mf,1)=t(nc,6)
          buff(2,m,mf,1)=t(nc,7)
!
       enddo
       if (bcg_to_proc(me)/=bcg_to_proc(other)) then
         call MPI_itrans2(buff(:,1:mt,mf,1),bcg_to_proc(me),bcg_to_proc(other),req(mf,1),me) ! send
         call MPI_itrans2(buff(:,1:mt,mf,2),bcg_to_proc(other),bcg_to_proc(me),req(mf,2),other) ! recv
       endif
    enddo

    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
       other=ndcc(mfb)
       me=bcl_to_bcg(mfb)
       if (bcg_to_proc(me)/=bcg_to_proc(other)) then
         call WAIT_MPI(req(mf,2))  ! waiting for the message to be received
       else
         buff(:,1:mt,mf,2)=buff(:,1:mt,bcg_to_mf(other),1)
       endif
!
       do m=1,mt
          mb=mpb(mfb)+m
          nd=ncbd(mb)
          ndm=ncin(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
          t(nd,6) = 0.5*( t(ndm,6)+buff(1,m,mf,2) )
          t(nd,7) = 0.5*( t(ndm,7)+buff(2,m,mf,2) )
!
       enddo
!
    enddo

    do mf=1,nbd
       mfb=lbd(mf)
       other=ndcc(mfb)
       me=bcl_to_bcg(mfb)
       if (bcg_to_proc(me)/=bcg_to_proc(other)) &
           call WAIT_MPI(req(mf,1))  ! waiting for all the messages to be sent
    enddo
    deallocate(buff)
!
    return
  end subroutine met_rbvc
end module mod_met_rbvc
