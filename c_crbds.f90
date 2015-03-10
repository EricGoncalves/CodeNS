      subroutine c_crbds( &
                 mot,imot,nmot, &
                 ncbd)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action crbds.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use sortiefichier
implicit none
integer :: imot
integer :: nmot
integer :: ncbd
integer :: imax
integer :: imin
integer :: jmax
integer :: jmin
integer :: kini
integer :: kmax
integer :: kmin
integer :: l
integer :: mfbe
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
      character(len=2 ) :: indmf
!
      dimension imot(nmx)
      dimension ncbd(ip41)
!
      call tcmd_crbds( &
                 mot,imot,nmot, &
                 mfbe,kini,l, &
                 imin,imax,jmin,jmax,kmin,kmax, &
                 indmf)
!
      if (kimp.ge.1) then
        call b1_crbds( &
                 mfbe,kini,l,imin,imax,jmin,jmax,kmin,kmax, &
                 indmf)
      endif
!
      call crbds( &
                 mfbe,kini,l, &
                 imin,imax,jmin,jmax,kmin,kmax, &
                 indmf, &
                 ncbd)
!
      if(kimp.ge.2) then
        call b2_crbds(mfbe)
      endif
!
      return
      end
