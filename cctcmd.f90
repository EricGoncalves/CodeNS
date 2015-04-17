module mod_cctcmd
  implicit none
contains
  subroutine cctcmd(command,lgcmd,mot,imot,nm1,nm2)
!
!***********************************************************************
!
!     ACT
!_A    Transformation d'une suite de mots en
!_A    une chaine de caractere.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    implicit none
    integer          ::      icmd,imot(nmx),     ipos,    lgcmd,       nm
    integer          ::       nm1,      nm2
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: command
    character(len=32) ::  mot(nmx)
!
    do icmd=1,lgcmdx
       command(icmd:icmd)=' '
    enddo
!
    lgcmd=0
    do nm=nm1,nm2
       lgcmd=lgcmd+1
       command(lgcmd:lgcmd)=' '
!$OMP SIMD
       do ipos=1,imot(nm)
          lgcmd=lgcmd+1
          command(lgcmd:lgcmd)=mot(nm)(ipos:ipos)
       enddo
    enddo
!
    return
  end subroutine cctcmd
end module mod_cctcmd
