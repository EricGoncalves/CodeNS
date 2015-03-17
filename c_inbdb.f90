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
    use para_fige
    use sortiefichier
    use mod_tcmd_inbdb
    use mod_b1_inbdb
    use mod_inbdb
    implicit none
    integer          :: ibdcfl,ibdcst,ibddim,  imot, kibdb
    integer          ::      l,  lmfb, lmfbd,  ncbd,  ncin
    integer          ::   nmot,  nvbc
    double precision :: bceqt,  vbc
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    character(len=4 ) :: clmf
!
    dimension imot(nmx)
    dimension lmfb(mtb)
    dimension vbc(ista*lsta)
!
    dimension ncbd(ip41)
    dimension ncin(ip41)
    dimension bceqt(ip41,neqt)
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
       call inbdb( &
            ncbd,ncin, &
            lmfb(l),clmf,kibdb, &
            ibdcst,ibdcfl,ibddim,nvbc,vbc,bceqt)
!
    enddo
!
    return
  end subroutine c_inbdb
end module mod_c_inbdb
