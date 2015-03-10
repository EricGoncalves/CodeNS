      subroutine c_crdms(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action crdms.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use sortiefichier
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
!
      call tcmd_crdms( &
                 mot,imot,nmot, &
                 l,ni,nj,nk)
!
      if (kimp.ge.1) then
            call b1_crdms(l,ni,nj,nk)
      endif
!
      call crdms(l,ni,nj,nk)
!
      if(kimp.ge.2) then
         call b2_crdms(l)
      endif
      if(kimp.ge.1) then
         call b3_crdms(l,ni,nj,nk)
      endif
!
      return
      end
