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
    select case(kvar)
    case(0)
       cvar=c0
    case(1)
       cvar=c1
    case(2)
       cvar=c2
    case(3)
       cvar=c3
    end select
!
    return
  end subroutine convich
end module mod_convich
