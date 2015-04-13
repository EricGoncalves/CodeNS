module mod_defdfpmdtd
  implicit none
contains
  subroutine defdfpmdtd
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfpmdtd.
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
       eta(l) =1.0

       keta(l)=1
    enddo
!
!$OMP END MASTER
    return
  end subroutine defdfpmdtd
end module mod_defdfpmdtd
