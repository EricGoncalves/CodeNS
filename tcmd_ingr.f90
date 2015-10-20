module mod_tcmd_ingr
  implicit none
contains
  subroutine tcmd_ingr( &
       mot,imot,nmot, &
       ldom,ldomd,king)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action ingr.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use maillage
    use kcle
    use mod_valenti
    use mod_vallent
    use mod_mpi
    use sortiefichier
    implicit none
    integer          ::       icmt, imot(nmx),      king,      kval,ldom(nobj)
    integer          ::      ldomd,        nm,      nmot,lzx2
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
       call sum_mpi(lzx,lzx2)
       write(stderr,*) rank,lzx2,klzx
       call vallent(mot,imot,nm,ldom,ldomd,lzx2,klzx)
    endif
!
    nm=nm+1
    if(nmot.lt.nm) then
       comment=ci
       call synterr(mot,imot,nmot,comment)
    else
       call valenti(mot,imot,nm,king,kval)
    endif
!
    return
  end subroutine tcmd_ingr
end module mod_tcmd_ingr
