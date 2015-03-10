module mod_defdfpmcfg
implicit none
contains
      subroutine defdfpmcfg
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfpmcfg.
!
!     COM
!_C    Par defaut les frontieres sont traites dans l'ordre dans lequel
!_C    elles ont ete creees.
!
!***********************************************************************
!
      use para_fige
   use kcle
   use maillage
   use boundary
implicit none
integer :: m
!
!-----------------------------------------------------------------------
!
      do m=1,mtb
       nba(m)=m
      enddo
!
      knba=1
!
      return
      end
end module
