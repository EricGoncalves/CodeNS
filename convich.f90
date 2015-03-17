module mod_convich
  implicit none
contains
  subroutine convich(kvar,cvar)
!
!***********************************************************************
!
!     ACT
!_A    Traduction d'une cle entiere(i) en
!_A    un commentaire (ci).
!
!***********************************************************************
!
    use chainecarac 
    implicit none
    integer          :: kvar
    character(len=24) ::  cvar
!
    if(kvar.eq.0) then
       cvar=c0
    else if (kvar.eq.1) then
       cvar=c1
    else if (kvar.eq.2) then
       cvar=c2
    else if (kvar.eq.3) then
       cvar=c3
    endif
!
    return
  end subroutine convich
end module mod_convich
