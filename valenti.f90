module mod_valenti
  use mod_synterr
  implicit none
contains
  subroutine valenti(mot,imot,nm,ient,kient)
!
!***********************************************************************
!
!     ACT
!_A    Affectation de sa valeur a l'entier ient.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use chainecarac
    use mod_entier
    implicit none
    integer          ::      icmt,     ient,imot(nmx),     kerr,    kient
    integer          ::        nm
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
!
!$OMP SIMD
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
!
    call entier(mot(nm),imot(nm),ient,kerr)
    if(kerr.eq.0)then
       comment=ci
       call synterr(mot,imot,nm,comment)
    else
       kient=2
    endif
!
    return
  end subroutine valenti
end module mod_valenti
