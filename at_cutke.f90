module mod_at_cutke
  implicit none
contains
  subroutine at_cutke(l,v)
!
!***********************************************************************
!
!     ACT
!_A
!_A   application des limiteurs aux variables turbulentes suivant
!_A   le modele choisi
!
!***********************************************************************
!
    use para_var
    use modeleturb
    use sortiefichier
    use mod_met_cutke
    use mod_met_cut
    use mod_met_cutsa
    use mod_met_cutked
    implicit none
    integer          :: l
    double precision :: v(ip11,ip60)
!
!-----------------------------------------------------------------------
!
!
    select case(kcutke)
    case(0)
       call met_cut(l,v)
    case(1)
!     les limiteurs sur k et epsilon lies ensembles
       call met_cutke(l,v)
    case(2)
!      Modele de Spalart Allmaras. Une seule equation (v[n,6])
       call met_cutsa(l,v)
    case(3)
!     les limiteurs sur k et epsilon sont decouples (famille k-omega)
       call met_cutked(l,v)
    case default
       write(imp,'(/,''!!!met_num: kcutke non prevu'')')
       stop
    end select
!
    return
  end subroutine at_cutke
end module mod_at_cutke
