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
implicit none
integer :: l
double precision :: v
!
!-----------------------------------------------------------------------
!
      dimension v(ip11,ip60)
!
      if(kcutke.eq.0) then
        call met_cut(l,v)
      elseif(kcutke.eq.1) then
!     les limiteurs sur k et epsilon lies ensembles
        call met_cutke(l,v)
      else if(kcutke.eq.2) then
!      Modele de Spalart Allmaras. Une seule equation (v[n,6])
        call met_cutsa(l,v)
      else if(kcutke.eq.3) then
!     les limiteurs sur k et epsilon sont decouples (famille k-omega)
        call met_cutked(l,v)
      else
        write(imp,'(/,''!!!met_num: kcutke non prevu'')')
        stop
      endif
!
       return
       end
