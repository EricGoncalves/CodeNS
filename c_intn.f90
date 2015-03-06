      subroutine c_intn(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action intn.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
	  use sortiefichier
!
!-----------------------------------------------------------------------
!
      character *32 mot(nmx)
!
      dimension imot(nmx)
!
      call tcmd_intn(mot,imot,nmot)
!
      if(kimp.ge.1) then
            call b1_intn
      endif
!
      return
      end
