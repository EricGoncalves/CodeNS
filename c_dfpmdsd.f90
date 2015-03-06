      subroutine c_dfpmdsd(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfpmdsd.
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
      dimension imot(nmx)
      dimension ldom(nobj)
      dimension lgr(nobj)
!
      call tcmd_dfpmdsd( &
                 mot,imot,nmot, &
                 ldom,ldomd, &
                 lgr,lgrd)
!
      if(kimp.ge.1) then
            call b1_dfpmdsd( &
                 ldom,ldomd, &
                 lgr,lgrd)
      endif
!
      return
      end
