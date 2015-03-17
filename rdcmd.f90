module mod_rdcmd
  implicit none
contains
  subroutine rdcmd(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Lecture d'une commande et transformation en une suite de mots.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use mod_gtcmd
    use mod_splcmd
    implicit none
    integer          ::  imot,lgcmd, nmot
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: command
    character(len=32) ::  mot(nmx)
    dimension imot(nmx)
!
    call gtcmd(command,lgcmd)
    call splcmd(command,lgcmd,mot,imot,nmot)
!
    return
  end subroutine rdcmd
end module mod_rdcmd
