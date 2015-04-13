module mod_c_end
  implicit none
contains
  subroutine c_end(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action end.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use mod_eend

    use mod_b1_end

    implicit none
    integer          :: imot(nmx),     nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
!$OMP MASTER
!
    call b1_end
    call eend
!
!$OMP END MASTER
    return
  end subroutine c_end
end module mod_c_end
