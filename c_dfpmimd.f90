module mod_c_dfpmimd
  implicit none
contains
  subroutine c_dfpmimd(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfpmimd.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use mod_b1_dfpmimd
    use mod_tcmd_dfpmimd
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
    call tcmd_dfpmimd( &
         mot,imot,nmot, &
         ldom,ldomd, &
         lgr,lgrd)
!
    if(kimp.ge.1) then
       if (rank==0) then
       call b1_dfpmimd( &
            ldom,ldomd, &
            lgr,lgrd)
       endif
       call barrier
    endif
!
    deallocate(ldom,lgr)
    return
  end subroutine c_dfpmimd
end module mod_c_dfpmimd
