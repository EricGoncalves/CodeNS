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
!
!-----------------------------------------------------------------------
!
      character *32 mot(nmx)
      character *4 clmf
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
      end
