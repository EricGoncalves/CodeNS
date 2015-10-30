module mod_met_rbsc
  implicit none
contains
  subroutine met_rbsc( &
       s1x,s1y,s1z,s2x,s2y,s2z, &
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
    double precision :: s1x(ip12),s1y(ip12),s1z(ip12),s2x(ip12),s2y(ip12)
    double precision :: s2z(ip12)
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
    allocate(buff(6,mt,nbd,2))

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
          buff(1,m,mf,1)=s1x(nc)
          buff(2,m,mf,1)=s1y(nc)
          buff(3,m,mf,1)=s1z(nc)
          buff(4,m,mf,1)=s2x(nc)
          buff(5,m,mf,1)=s2y(nc)
          buff(6,m,mf,1)=s2z(nc)
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
          s1x(nd) = 0.5*( s1x(ndm)+buff(1,m,mf,2) )
          s1y(nd) = 0.5*( s1y(ndm)+buff(2,m,mf,2) )
          s1z(nd) = 0.5*( s1z(ndm)+buff(3,m,mf,2) )
          s2x(nd) = 0.5*( s2x(ndm)+buff(4,m,mf,2) )
          s2y(nd) = 0.5*( s2y(ndm)+buff(5,m,mf,2) )
          s2z(nd) = 0.5*( s2z(ndm)+buff(6,m,mf,2) )
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
  end subroutine met_rbsc
end module mod_met_rbsc
