module mod_c_crdms
  implicit none
contains
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
    use mod_crdms
    use mod_b3_crdms
    use mod_b1_crdms
    use mod_b2_crdms
    use mod_tcmd_crdms
    use mod_mpi,only : rank
    implicit none
    integer          :: imot(nmx),        l,       ni,       nj,       nk
    integer          ::      nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
!
    call tcmd_crdms( &
         mot,imot,nmot, &
         l,ni,nj,nk)

    if(l==rank+1) then
  !
      if (kimp.ge.1) then
         call b1_crdms(l,ni,nj,nk)
      endif
    endif
  !
      call crdms(l,ni,nj,nk)
  !
    if(l==rank+1) then
      if(kimp.ge.2) then
         call b2_crdms(l)
      endif
      if(kimp.ge.1) then
         call b3_crdms(l,ni,nj,nk)
      endif
    endif
!
    return
  end subroutine c_crdms
end module mod_c_crdms
