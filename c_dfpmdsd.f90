module mod_c_dfpmdsd
  implicit none
contains
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
    use mod_b1_dfpmdsd
    use mod_tcmd_dfpmdsd
    use mod_mpi
    implicit none
    integer          :: imot(nmx),    ldomd,     lgrd,     nmot
    integer         ,allocatable :: ldom(:), lgr(:)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    allocate(ldom(nobj))
    allocate(lgr(nobj))
!
    call tcmd_dfpmdsd( &
         mot,imot,nmot, &
         ldom,ldomd, &
         lgr,lgrd)
!
    if(kimp.ge.1) then
        if (rank==0) then
        call b1_dfpmdsd( &
            ldom,ldomd, &
            lgr,lgrd)
        endif
        call barrier
    endif
!
    deallocate(ldom,lgr)
    return
  end subroutine c_dfpmdsd
end module mod_c_dfpmdsd
