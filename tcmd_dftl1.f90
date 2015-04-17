module mod_tcmd_dftl1
  implicit none
contains
  subroutine tcmd_dftl1(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dftl1.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use kcle
    use mod_cctcmd
    implicit none
    integer          :: imot(nmx),     itit,    lgcmd,     nmot
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: command
    character(len=32)   :: mot(nmx)
!
    if(ktitrt1.eq.2) ktitrt1=3
!
    if(nmot.eq.2)then
       do itit=1,80
          titrt1(itit:itit)=' '
       enddo
    else
       call cctcmd(command,lgcmd,mot,imot,3,nmot)
!$OMP SIMD
       do itit=1,80
          titrt1(itit:itit)=command(itit:itit)
       enddo
    endif
!
    ktitrt1=2
!
    return
  end subroutine tcmd_dftl1
end module mod_tcmd_dftl1
