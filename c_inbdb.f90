module mod_c_inbdb
  implicit none
contains
  subroutine c_inbdb( &
       mot,imot,nmot, &
       ncbd,ncin,bceqt,partition)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action inbdb.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use boundary
    use sortiefichier
    use mod_tcmd_inbdb
    use mod_b1_inbdb
    use mod_inbdb
    use mod_mpi
    implicit none
    integer          ::     ibdcfl,    ibdcst,    ibddim, imot(nmx),     kibdb
    integer          ::          l,     lmfbd,ncbd(ip41),ncin(ip41),      nmot
    integer          ::       nvbc,l1
    double precision :: bceqt(ip41,neqt),  vbc(ista*lsta)
    integer         ,allocatable :: lmfb(:)
    logical,optional :: partition
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    character(len=4 ) :: clmf
    allocate(lmfb(mtb))
!
!
!
    call tcmd_inbdb( &
         mot,imot,nmot, &
         lmfb,lmfbd,clmf,kibdb, &
         ibdcst,ibdcfl,ibddim,nvbc,vbc)
!
    if (kimp.ge.1) then
       if (rank==0) call b1_inbdb( &
            lmfb,lmfbd,clmf,kibdb, &
            ibdcst,ibdcfl,ibddim,nvbc,vbc)
    endif
!
    do l=1,lmfbd
!
      if (present(partition)) then
       call inbdb( &
            ncbd,ncin, &
            lmfb(l),clmf,kibdb, &
            ibdcst,ibdcfl,ibddim,nvbc,vbc,bceqt)
    else
       do l1=1,mtb
        if (bcg_to_bci(bcl_to_bcg(l1))==lmfb(l)) &
       call inbdb( &
            ncbd,ncin, &
            l1,clmf,kibdb, &
            ibdcst,ibdcfl,ibddim,nvbc,vbc,bceqt)
!
    enddo
    endif
    enddo
!
    deallocate(lmfb)
    return
  end subroutine c_inbdb
end module mod_c_inbdb
