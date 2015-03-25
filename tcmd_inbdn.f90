module mod_tcmd_inbdn
  implicit none
contains
  subroutine tcmd_inbdn( &
       mot,imot,nmot, &
       lmfb,lmfbd,kibdn)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action inbdn.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use maillage
    use kcle
    use mod_valenti
    use mod_vallent
    implicit none
  integer          ::       icmt, imot(nmx),     kibdn,      kval,lmfb(nobj)
  integer          ::      lmfbd,        nm,      nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
!
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
    kval=0
!
    nm=3
    nm=nm+1
    if(nmot.lt.nm) then
       comment=cm
       call synterr(mot,imot,nmot,comment)
    else
       call vallent(mot,imot,nm,lmfb,lmfbd,mtbx,kmtbx)
    endif
!
    nm=nm+1
    if(nmot.lt.nm) then
       comment=ci
       call synterr(mot,imot,nmot,comment)
    else
       call valenti(mot,imot,nm,kibdn,kval)
    endif
!
    return
  end subroutine tcmd_inbdn
end module mod_tcmd_inbdn
