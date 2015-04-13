module mod_b1_dftl1
  implicit none
contains
  subroutine b1_dftl1
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dftl1'.
!
!     INP
!_I    titrt1     : com char             ; titre du calcul
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!
    use sortiefichier
    use chainecarac
    use kcle
    use mod_convich
    implicit none
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=24) ::  ctitrt1
!$OMP MASTER
!
    call convich(ktitrt1,ctitrt1)
!
    form='(/,2x,''titre du calcul : '',/' &
         //'2x,''---------------'',/' &
         //'2x,a,2x,a)'
    write(imp,form) titrt1,ctitrt1
!
!$OMP END MASTER
    return
  end subroutine b1_dftl1
end module mod_b1_dftl1
