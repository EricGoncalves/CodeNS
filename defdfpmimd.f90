module mod_defdfpmimd
  implicit none
contains
  subroutine defdfpmimd
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfpmimd.
!
!     COM
!_C    Par defaut utilisation de la phase implicite avec des parametres
!_C    imposes par defaut.
!
!***********************************************************************
!
    use para_fige
    use schemanum
    use kcle
    implicit none
    integer          :: l
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!$OMP SIMD
    do l=1,lz
       kmf(l) =1
       lmax(l)=10

       kkmf(l) =1
       klmax(l)=1
    enddo
!
!$OMP END MASTER
    return
  end subroutine defdfpmimd
end module mod_defdfpmimd
