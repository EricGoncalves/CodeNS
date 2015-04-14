module mod_tcmd_intn
  implicit none
contains
  subroutine tcmd_intn(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action intn.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use kcle
    use schemanum
    use mod_valenti
    implicit none
    integer          ::      icmt,imot(nmx),       nm,     nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
!
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
!
    if(knumt.eq.2) knumt=3
!
    if(nmot.eq.3)then
       comment=cb
       call synterr(mot,imot,2,comment)
    endif
!
    nm=3
!
    nm=nm+1
    if(nmot.lt.nm) then
       comment=ci
       call synterr(mot,imot,nmot,comment)
    else
       call valenti(mot,imot,nm,numt,knumt)
    endif
!
    return
  end subroutine tcmd_intn
end module mod_tcmd_intn
