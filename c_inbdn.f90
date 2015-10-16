module mod_c_inbdn
  implicit none
contains
  subroutine c_inbdn( &
       mot,imot,nmot, &
       x,y,z, &
       sn, &
       ncbd,nxn,nyn,nzn, &
       tn1,tn2,tn3,tn4,tn5,tn6, &
       tn7,tn8,tn9)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action inbdn.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use boundary
    use sortiefichier
    use mod_b1_inbdn
    use mod_inbdn
    use mod_tcmd_inbdn
    use mod_mpi
    implicit none
    integer          ::  imot(nmx),     kibdn,         l,     lmfbd,ncbd(ip41)
    integer          ::       nmot,l1,l2
    double precision ::     nxn(ip42),    nyn(ip42),    nzn(ip42),sn(ip31*ndir),    tn1(ip00)
    double precision ::     tn2(ip00),    tn3(ip00),    tn4(ip00),    tn5(ip00),    tn6(ip00)
    double precision ::     tn7(ip00),    tn8(ip00),    tn9(ip00),      x(ip21),      y(ip21)
    double precision ::       z(ip21)
    integer         ,allocatable :: lmfb(:)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    allocate(lmfb(mtb))
!
!
!
!
!
    call tcmd_inbdn( &
         mot,imot,nmot, &
         lmfb,lmfbd,kibdn)
!
    if (kimp.ge.1) then
       if (rank==0) call b1_inbdn(lmfb,lmfbd,kibdn)
    endif
!
    do l=1,lmfbd
!
       do l1=1,mtb
        l2=bcl_to_bcg(l1)
        if (bcg_to_bci(l2)==lmfb(l)) &
       call inbdn( &
            l2,kibdn, &
            x,y,z, &
            sn, &
            ncbd,nxn,nyn,nzn, &
            tn1,tn2,tn3,tn4,tn5,tn6, &
            tn7,tn8,tn9)
   enddo
!
    enddo
!
    deallocate(lmfb)
    return
  end subroutine c_inbdn
end module mod_c_inbdn
