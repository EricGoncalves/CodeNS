      subroutine c_secpfw(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action secpfw.
!
!***********************************************************************
!
      use para_fige
   use sortiefichier
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
!
      dimension imot(nmx)
      dimension lgr(nobj)
!
      call tcmd_secpfw( &
                 mot,imot,nmot, &
                 lgr,lgrd)
!
      if(kimp.ge.1) then
            call b1_secpfw(lgr,lgrd)
      endif
!
      return
      end
