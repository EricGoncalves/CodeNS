module mod_valreel
  use mod_synterr
  implicit none
contains
  subroutine valreel(mot,imot,nm,rree,krree)
!
!***********************************************************************
!
!     ACT
!_A    Affectation de sa valeur au reel rree.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use mod_reel
    implicit none
    integer          ::      icmt,imot(nmx),     kerr,    krree,       nm
    double precision :: rree
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
!$OMP MASTER
!
!$OMP SIMD
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
!
    call reel(mot(nm),imot(nm),rree,kerr)
    if(kerr.eq.0)then
       comment=cr
       call synterr(mot,imot,nm,comment)
    else
       krree=2
    endif
!
!$OMP END MASTER
    return
  end subroutine valreel
end module mod_valreel
