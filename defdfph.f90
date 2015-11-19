module mod_defdfph
  implicit none
contains
  subroutine defdfph
!
!***********************************************************************
!
!_DA  DATE_C : mai 2006 -- Eric GONCALVES / LEG
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfph.
!
!      pour air
!_C    Par defaut gamma=1.4 et constante des gaz parfait r=287.13
!_C    Par defaut Prandtl de 0,72 et Prandtl turbulent de 0,9 .
!
!      pour eau
!_C    Par defaut gamma=1.01 et constante des gaz parfait r=41,4158
!_C    Par defaut Prandtl de 7 et Prandtl turbulent de 1.
!
!***********************************************************************
!
    use kcle
    use constantes
    use proprieteflu
    implicit none
!
!-----------------------------------------------------------------------
!
    gam=1.4 
    kgam=1
!
    rd=287.13
!    rd=41.4158
    krd=1
!
    pr=0.72
    kpr=1
!
    prt=1.
    kprt=1
!
    reynz=reelmx
    kreynz=0
!
    return
  end subroutine defdfph
end module mod_defdfph
