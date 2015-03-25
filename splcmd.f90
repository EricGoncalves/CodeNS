module mod_splcmd
  implicit none
contains
  subroutine splcmd(command,lgcmd,mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Transformation d'une chaine de caracteres en
!_A    une suite de mots par detection des blancs.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    implicit none
  integer          ::      icmd,       im,imot(nmx),      kbl,    lgcmd
  integer          ::        nm,     nmot
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: command
    character(len=32) ::  mot(nmx)
!
!
    nmot=0
    do nm=1,50
       do im=1,32
          mot(nm)(im:im)=' '
          imot(nm)=0
       enddo
    enddo
    kbl=0
    im=0
!
    do icmd=1,lgcmd
       if (command(icmd:icmd).eq.' ') then
          if((nmot.ge.1).and.(kbl.eq.0)) then
             imot(nmot)=im
             nmot=nmot+1
             im=0
             kbl=1
          endif
       else
          if(nmot.eq.0) nmot=1
          kbl=0
          im=im+1
          mot(nmot)(im:im)=command(icmd:icmd)
       endif
    enddo
    if (nmot.ne.0) imot(nmot)=im
!
    return
  end subroutine splcmd
end module mod_splcmd
