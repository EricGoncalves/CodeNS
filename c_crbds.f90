module mod_c_crbds
  implicit none
contains
  subroutine c_crbds( &
       mot,imot,nmot, &
       ncbd,mfbi)
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
    use mod_b2_crbds
    use mod_b1_crbds
    use mod_tcmd_crbds
    use mod_crbds
    use mod_mpi
    implicit none
    integer          ::       imax,      imin, imot(nmx),      jmax,      jmin
    integer          ::       kini,      kmax,      kmin,         l,      mfbe
    integer          ::       nmot
    integer,optional :: mfbi
    integer,allocatable :: ncbd(:)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    character(len=2 ) :: indmf
!
!
    call tcmd_crbds( &
         mot,imot,nmot, &
         mfbe,kini,l, &
         imin,imax,jmin,jmax,kmin,kmax, &
         indmf)
!
    if (kimp.ge.1) then
      if(rank==0) then
         call b1_crbds( &
              mfbe,kini,l,imin,imax,jmin,jmax,kmin,kmax, &
              indmf)
      endif
      call barrier
    endif
!
    if (present(mfbi)) then
      call crbds( &
           mfbe,kini,l, &
           imin,imax,jmin,jmax,kmin,kmax, &
           indmf, &
           ncbd,mfbi)
    else
      call crbds( &
           mfbe,kini,l, &
           imin,imax,jmin,jmax,kmin,kmax, &
           indmf, &
           ncbd)
    endif

!
    if(bcg_to_proc(mfbe)==rank) then
      if(kimp.ge.2) then
         call b2_crbds(bcg_to_bcl(mfbe))
      endif
    endif
!
    return
  end subroutine c_crbds
end module mod_c_crbds
