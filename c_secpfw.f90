module mod_c_secpfw
  implicit none
contains
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
    use mod_tcmd_secpfw
    use mod_b1_secpfw
    use mod_mpi
    implicit none
    integer          :: imot(nmx),     lgrd,     nmot
    integer         ,allocatable :: lgr(:)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    allocate(lgr(nobj))
!
!
    call tcmd_secpfw( &
         mot,imot,nmot, &
         lgr,lgrd)
!
    if(kimp.ge.1) then
       if (rank==0) call b1_secpfw(lgr,lgrd)
       call barrier
    endif
!
    deallocate(lgr)
    return
  end subroutine c_secpfw
end module mod_c_secpfw
