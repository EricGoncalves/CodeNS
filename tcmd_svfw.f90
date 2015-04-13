module mod_tcmd_svfw
  implicit none
contains
  subroutine tcmd_svfw( &
       mot,imot,nmot, &
       disc)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action svfw.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use mod_synterr
    implicit none
    integer          ::      icmt,       im,imot(nmx),       nm,     nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
    character(len=4 ) :: disc
!$OMP MASTER
!
!$OMP SIMD
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
!
    nm=3
    nm=nm+1
    if(nmot.lt.nm) then
       comment=ch
       call synterr(mot,imot,nmot,comment)
    else
       if(imot(nm).gt.4)then
          comment=cc
          call synterr(mot,imot,nm,comment)
       else
          disc(1:4)='    '
!$OMP SIMD
          do im=1,imot(nm)
             disc(im:im)=mot(nm)(im:im)
          enddo
       endif
    endif
!
!$OMP END MASTER
    return
  end subroutine tcmd_svfw
end module mod_tcmd_svfw
