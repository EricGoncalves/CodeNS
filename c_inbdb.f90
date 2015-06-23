module mod_c_inbdb
  implicit none
contains
  subroutine c_inbdb( &
       mot,imot,nmot, &
       ncbd,ncin,bceqt)
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
    implicit none
    integer          ::     ibdcfl,    ibdcst,    ibddim, imot(nmx),     kibdb
    integer          ::          l,     lmfbd,ncbd(ip41),ncin(ip41),      nmot
    integer          ::       nvbc,l1
    double precision :: bceqt(ip41,neqt),  vbc(ista*lsta)
    integer         ,allocatable :: lmfb(:)
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
       call b1_inbdb( &
            lmfb,lmfbd,clmf,kibdb, &
            ibdcst,ibdcfl,ibddim,nvbc,vbc)
    endif
!
    do l=1,lmfbd
!
       do l1=1,mtb
        if (new2old_f(l1)==lmfb(l)) &
       call inbdb( &
            ncbd,ncin, &
            l1,clmf,kibdb, &
            ibdcst,ibdcfl,ibddim,nvbc,vbc,bceqt)
!
    enddo
    enddo
!
    deallocate(lmfb)
    return
  end subroutine c_inbdb
end module mod_c_inbdb
