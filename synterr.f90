module mod_synterr
  implicit none
contains
  subroutine synterr(mot,imot,nmot,comment)
!
!***********************************************************************
!
!     ACT
!_A    Message d'erreur en cas de donnees.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use mod_cctcmd
    implicit none
    integer          :: imot(nmx),     ipos,    lgcmd,     nmot
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: command
    character(len=32) :: comment
    character(len=32) :: mot(nmx)
    character(len=7 ) :: formatcm
    character(len=4 ) :: longcm
!$OMP MASTER
!
    call cctcmd(command,lgcmd,mot,imot,1,nmot)
!
    if (lgcmd.lt.(lgcmdx-3)) then
       lgcmd=lgcmd+1
       command(lgcmd:lgcmd)='   '
       do ipos=1,3
          lgcmd=lgcmd+1
          command(lgcmd:lgcmd)='.'
       enddo
    endif
!
    write(longcm,'(i4)') lgcmd
    formatcm='(a'//longcm//')'
!
    write(imp,'(/a)') &
         ' !! Attention !! erreur detectee dans la commande :'
    write(imp,formatcm) command
    write(imp,'(a)') comment
!
    stop 'Erreur de syntaxe dans une commande!'
!
!$OMP END MASTER
    return
  end subroutine synterr
end module mod_synterr
